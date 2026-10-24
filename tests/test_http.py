import json
import tempfile
import threading
import unittest
import urllib.error
import urllib.request
from http.server import ThreadingHTTPServer
from pathlib import Path

from spatial_lab import server
from spatial_lab.data import digest


class HttpTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.old = server.RUNTIME
        server.RUNTIME = Path(cls.tmp.name)
        cls.http = ThreadingHTTPServer(("127.0.0.1", 0), server.Handler)
        cls.thread = threading.Thread(target=cls.http.serve_forever, daemon=True)
        cls.thread.start()
        cls.base = f"http://127.0.0.1:{cls.http.server_port}"

    @classmethod
    def tearDownClass(cls):
        cls.http.shutdown()
        cls.http.server_close()
        cls.thread.join()
        server.RUNTIME = cls.old
        cls.tmp.cleanup()

    def request(self, path, body=None, token=True, origin=None):
        headers = {}
        if body is not None:
            headers["Content-Type"] = "application/json"
            if token:
                headers["X-Lab-Token"] = server.TOKEN
        if origin:
            headers["Origin"] = origin
        r = urllib.request.Request(
            self.base + path,
            headers=headers,
            data=json.dumps(body).encode() if body is not None else None,
        )
        with urllib.request.urlopen(r, timeout=10) as response:
            return json.load(response)

    def test_bootstrap_and_region_have_provenance(self):
        b = self.request("/api/bootstrap")
        r = self.request("/api/region?region=reg009_A")
        self.assertEqual(b["summary"]["rows"], 258385)
        self.assertEqual(r["metadata"]["patient"], 5)
        self.assertEqual(len(r["cell_id"]), 2434)
        self.assertIn("Not assigned by authors", r["neighborhood_names"])

    def test_reconstruction_record_matches_export(self):
        r = self.request(
            "/api/run",
            {
                "kind": "reconstruction",
                "region": "reg009_A",
                "method": "ridge",
                "pattern": "band",
                "seed": 107,
            },
        )
        self.assertEqual(digest(r["result"]), r["result_sha256"])
        self.assertEqual(r["parameters"]["pattern"], "band")
        self.assertEqual(r, self.request("/api/run/" + r["id"]))
        d = self.request("/api/dossier?region=reg009_A")
        self.assertIn(r, d["history"])
        self.assertEqual(d["data"]["derived_sha256"], r["data_sha256"])
        self.assertNotIn("token", d)

    def test_rejections_do_not_create_history(self):
        before = self.request("/api/history")
        for body in [
            {"kind": "movement", "region": "reg009_A", "permeability": True},
            {"kind": "movement", "region": "reg009_A", "seed": True},
            {"kind": "reconstruction", "region": "unknown"},
            {"kind": "note", "region": "reg009_A", "note": ""},
        ]:
            with self.assertRaises(urllib.error.HTTPError) as caught:
                self.request("/api/run", body)
            self.assertEqual(caught.exception.code, 400)
        self.assertEqual(before, self.request("/api/history"))

    def test_foreign_origin_and_missing_token_are_rejected(self):
        for kwargs in [{"token": False}, {"origin": "https://unrelated.example"}]:
            with self.assertRaises(urllib.error.HTTPError) as caught:
                self.request(
                    "/api/run", {"kind": "sensitivity", "region": "reg009_A"}, **kwargs
                )
            self.assertEqual(caught.exception.code, 403)
        with self.assertRaises(urllib.error.HTTPError) as caught:
            self.request("/../data/source/CRC_clusters_neighborhoods_markers.csv")
        self.assertEqual(caught.exception.code, 404)

    def test_note_and_simulation_are_preserved_separately(self):
        note = self.request(
            "/api/run",
            {
                "kind": "note",
                "region": "reg009_A",
                "note": "Test interpretation: proximity is not cytotoxicity.",
            },
        )
        motion = self.request(
            "/api/run",
            {
                "kind": "movement",
                "region": "reg009_A",
                "permeability": 0.1,
                "seed": 101,
            },
        )
        self.assertEqual(motion["result"]["agents"], 64)
        self.assertEqual(len(motion["result"]["paths"]), 64)
        self.assertNotEqual(note["id"], motion["id"])


if __name__ == "__main__":
    unittest.main()
