from http import server
import socketserver
import threading
from contextlib import contextmanager
from urllib.parse import unquote
from pathlib import Path

class FileHandler(server.SimpleHTTPRequestHandler):

    def __init__(self, *args, file_path: str|Path, hosted_path: str, **kwargs):

        self.file_path = Path(file_path)
        self.hosted_path = "/" + hosted_path.lstrip("/")
        super().__init__(*args, **kwargs)

    def do_GET(self):
        if unquote(self.path) != self.hosted_path:
            self.send_error(404, "File not found")

        else:
            self.send_response(200)
            self.send_header("Content-type", "application/octet-stream")
            self.send_header("Content-Length", str(self.file_path.stat().st_size))
            self.end_headers()

            with open(self.file_path, "rb") as f:
                self.wfile.write(f.read())



class LocalFile:

    def __init__(self, file_path: str|Path, port: int=8007):

        self.file_path = Path(file_path).resolve()
        if not self.file_path.exists():
            raise FileNotFoundError(self.file_path)

        self.port = port
        self.hosted_path = self.file_path.name
        self.url = f"http://localhost:{self.port}/{self.hosted_path}"
        self._httpd = None
        self._thread = None

    @contextmanager
    def serve(self):
        def handler(*args, **kwargs):
            return FileHandler(*args, file_path=self.file_path, hosted_path=self.hosted_path, **kwargs)

        with socketserver.TCPServer(("", self.port), handler) as httpd:
            self._httpd = httpd
            self._thread = threading.Thread(target=httpd.serve_forever, daemon=True)
            self._thread.start()

            try:
                yield self.url
            finally:
                httpd.shutdown()
                self._thread.join()


