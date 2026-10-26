import gzip
import io
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import threading
from pathlib import Path

import pytest

from src.stream_extractor import remote_operation, HashReader, bounded_lines, aggregate_fragments, Blacklist, RemoteResourceError


def test_real_http_multimember_bgzf_style_and_retry_without_double_count(reference_bins, extract_params, monkeypatch):
    # gzip concatenated members exercise the same decompression contract as BGZF.
    payload = gzip.compress(b'chr1\t0\t150\t60\t+\n') + gzip.compress(b'chr1\t0\t151\t60\t-\n')
    class Handler(BaseHTTPRequestHandler):
        attempts = 0
        def do_GET(self):
            Handler.attempts += 1
            self.send_response(503 if Handler.attempts == 1 else 200)
            self.end_headers()
            self.wfile.write(b'busy' if Handler.attempts == 1 else payload)
        def log_message(self, *args):
            pass
    server = ThreadingHTTPServer(('127.0.0.1', 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        def consume(response):
            reader = HashReader(response.raw)
            with gzip.GzipFile(fileobj=reader) as handle:
                result = aggregate_fragments(bounded_lines(handle), reference_bins, extract_params, Blacklist())
            return result, reader.count
        (values, qc), size = remote_operation(f'http://127.0.0.1:{server.server_port}/test.frag.tsv.bgz',
                     dict(max_retries=5, backoff_seconds=0, connect_timeout=2, read_timeout=2), consume)
        assert Handler.attempts == 2 and qc['accepted'] == 2 and size == len(payload)
        assert values['ratio__chr1:0-5000000'] == 1
    finally:
        server.shutdown()
        server.server_close()
        thread.join()


def test_remote_failure_has_url_and_status_and_exact_retry_limit():
    class Response:
        status_code = 404
        def __enter__(self): return self
        def __exit__(self, *args): pass
        def raise_for_status(self):
            import requests
            raise requests.HTTPError('missing')
    class Session:
        count = 0
        def get(self, *args, **kwargs):
            self.count += 1
            return Response()
    s = Session()
    with pytest.raises(RemoteResourceError, match='URL=https://example.org/missing HTTP=404'):
        remote_operation('https://example.org/missing', dict(max_retries=5, backoff_seconds=0,
                         connect_timeout=1, read_timeout=1), lambda _: None, session=s)
    assert s.count == 6
