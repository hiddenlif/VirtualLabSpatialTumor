"""Loopback-only local research interface. Scientific computations run serially."""

import argparse
import csv
import io
import json
import re
import secrets
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, urlparse

from .analysis import grid, sensitivity
from .data import MARKERS, ROOT, digest, load, region
from .movement import simulate
from .reconstruction import frozen_model, reconstruct

TOKEN = secrets.token_urlsafe(24)
RUNTIME = ROOT / ".runtime"
LOCK = threading.Lock()
COMPUTE = threading.BoundedSemaphore(1)


def history():
    if not RUNTIME.exists():
        return []
    records = []
    for p in RUNTIME.glob("*.json"):
        value = json.loads(p.read_text())
        record = {
            k: value[k]
            for k in (
                "id",
                "kind",
                "region",
                "created_at",
                "parameters",
                "result_sha256",
            )
        }
        record["created_ns"] = value.get("created_ns", 0)
        records.append(record)
    return sorted(
        records, key=lambda x: (x["created_at"], x["created_ns"]), reverse=True
    )


def retained():
    folder = ROOT / "reports/study-v1"
    return {
        k: json.loads((folder / (name + ".json")).read_text())
        for k, name in [
            ("summary", "summary"),
            ("regions", "regions"),
            ("sensitivity", "sensitivity"),
        ]
    }


def run(body):
    if not isinstance(body, dict):
        # Keep malformed bodies on the public ValueError / HTTP 400 validation path.
        raise ValueError("Expected an object")  # noqa: TRY004
    name = body.get("region")
    r = region(name)
    kind = body.get("kind")
    seed = body.get("seed", 101)
    if type(seed) is not int or not 0 <= seed <= 1000000:
        raise ValueError("Seed must be an integer in 0..1000000")
    if kind == "reconstruction":
        parameters = {
            "method": body.get("method", "ridge"),
            "pattern": body.get("pattern", "random"),
            "seed": seed,
        }
        result = reconstruct(r, **parameters, include_cells=True)
    elif kind == "movement":
        permeability = body.get("permeability", 0.5)
        if type(permeability) not in (float, int):
            raise ValueError("Invalid permeability")
        parameters = {"permeability": permeability, "seed": seed}
        result = simulate(r, **parameters)
    elif kind == "sensitivity":
        parameters = {"seed": seed}
        result = {"region": name, "scenarios": sensitivity(r, seed)}
    elif kind == "note":
        note = body.get("note")
        if not isinstance(note, str) or not 1 <= len(note.strip()) <= 2000:
            raise ValueError("Interpretation must contain 1..2000 characters")
        parameters = {}
        result = {
            "region": name,
            "note": note.strip(),
            "author": "Local user interpretation; not an automated scientific conclusion",
        }
    else:
        raise ValueError("Unknown operation")
    record = {
        "id": secrets.token_hex(6),
        "kind": kind,
        "region": name,
        "created_at": time.strftime("%Y-%m-%d %H:%M:%S %Z"),
        "created_ns": time.time_ns(),
        "parameters": parameters,
        "data_sha256": load()[1]["derived_sha256"],
        "model_sha256": digest(frozen_model()),
        "result_sha256": digest(result),
        "result": result,
    }
    RUNTIME.mkdir(exist_ok=True)
    with LOCK:
        target = RUNTIME / (record["id"] + ".json")
        if target.exists():
            raise FileExistsError("Record ID collision")
        pending = RUNTIME / ("." + record["id"] + ".pending")
        with pending.open("x") as f:
            json.dump(record, f, allow_nan=False)
        pending.replace(target)
    return record


def region_payload(name):
    r = region(name)
    m = load()[1]
    return {
        "region": name,
        "metadata": r["metadata"],
        "markers": MARKERS,
        "types": m["type_names"],
        "neighborhood_names": m["neighborhood_names"],
        "cell_id": r["cell_id"].tolist(),
        "xy": r["xy"].tolist(),
        "labels": r["labels"].tolist(),
        "eligible": r["eligible"].tolist(),
        "marker_values": r["markers"].tolist(),
        "neighborhood": r["neighborhood"].tolist(),
        "terrain": grid(r),
        "inference": "Author cell identities are inferred annotations. Tumor-rich territory is this project’s k=20 neighborhood inference. Points are processed cell centroids, not an image.",
    }


class Handler(BaseHTTPRequestHandler):
    def log_message(self, *args):
        pass

    def valid_host(self):
        port = self.server.server_port
        return self.headers.get("Host") in (f"127.0.0.1:{port}", f"localhost:{port}")

    def send_bytes(self, b, status=200, content_type="application/json", filename=None):
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(b)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header(
            "Content-Security-Policy",
            "default-src 'self'; script-src 'self'; style-src 'self' 'unsafe-inline'; img-src 'self' data:; connect-src 'self'; object-src 'none'; frame-ancestors 'none'",
        )
        if filename:
            self.send_header(
                "Content-Disposition", f'attachment; filename="{filename}"'
            )
        self.end_headers()
        self.wfile.write(b)

    def json(self, obj, status=200):
        self.send_bytes(json.dumps(obj, allow_nan=False).encode(), status)

    def do_GET(self):
        if not self.valid_host():
            return self.json({"error": "Host rejected"}, 403)
        u = urlparse(self.path)
        q = parse_qs(u.query)
        name = q.get("region", [""])[0]
        try:
            if u.path == "/api/bootstrap":
                _, m = load()
                study = retained()
                return self.json(
                    {
                        "token": TOKEN,
                        "manifest": m,
                        "summary": study["summary"],
                        "regions": study["regions"],
                        "history": history(),
                    }
                )
            if u.path == "/api/region":
                return self.json(region_payload(name))
            if u.path == "/api/study":
                return self.json(retained())
            if u.path == "/api/history":
                return self.json(history())
            if u.path == "/api/dossier":
                entries = [
                    json.loads((RUNTIME / (x["id"] + ".json")).read_text())
                    for x in history()
                    if x["region"] == name
                ]
                r = region(name)
                return self.json(
                    {
                        "schema_version": "1.0",
                        "region": name,
                        "patient": r["patient"],
                        "data": load()[1],
                        "reference": next(
                            x for x in retained()["regions"] if x["region"] == name
                        ),
                        "history": entries,
                        "limitations": [
                            "2D processed centroids; no raw-image segmentation validation",
                            "Marker reconstruction is evaluated against retained measured intensities",
                            "Movement is hypothetical; no clinical or therapeutic outcome prediction",
                        ],
                    }
                )
            if u.path == "/api/cells.csv":
                r = region(name)
                m = load()[1]
                out = io.StringIO()
                w = csv.writer(out)
                w.writerow(
                    [
                        "cell_id",
                        "region",
                        "patient",
                        "x_source_px",
                        "y_source_px",
                        "author_cell_type",
                        "eligible",
                        "CD8_intensity_au",
                        "author_neighborhood",
                    ]
                )
                for i, c in enumerate(r["cell_id"]):
                    w.writerow(
                        [
                            int(c),
                            name,
                            r["patient"],
                            *r["xy"][i],
                            r["labels"][i],
                            bool(r["eligible"][i]),
                            r["markers"][i, 0],
                            m["neighborhood_names"][r["neighborhood"][i]],
                        ]
                    )
                return self.send_bytes(
                    out.getvalue().encode(),
                    content_type="text/csv",
                    filename=name + "-cells.csv",
                )
            if u.path.startswith("/api/run/"):
                key = u.path.rsplit("/", 1)[-1]
                if not re.fullmatch("[0-9a-f]{12}", key):
                    raise ValueError("Unknown record")
                return self.json(json.loads((RUNTIME / (key + ".json")).read_text()))
            files = {
                "/": ("index.html", "text/html"),
                "/app.js": ("app.js", "application/javascript"),
                "/style.css": ("style.css", "text/css"),
            }
            if u.path not in files:
                return self.json({"error": "Not found"}, 404)
            path, ctype = files[u.path]
            return self.send_bytes(
                (ROOT / "web" / path).read_bytes(), content_type=ctype
            )
        except (ValueError, KeyError, FileNotFoundError, TypeError) as e:
            return self.json({"error": str(e)}, 400)

    def do_POST(self):
        if not self.valid_host():
            return self.json({"error": "Host rejected"}, 403)
        allowed = (
            f"http://127.0.0.1:{self.server.server_port}",
            f"http://localhost:{self.server.server_port}",
        )
        if (
            self.headers.get("X-Lab-Token") != TOKEN
            or self.headers.get("Origin", allowed[0]) not in allowed
        ):
            return self.json({"error": "Local origin and request token required"}, 403)
        if self.path != "/api/run":
            return self.json({"error": "Not found"}, 404)
        if not COMPUTE.acquire(blocking=False):
            return self.json({"error": "A computation is already running"}, 429)
        try:
            length = int(self.headers.get("Content-Length", "0"))
            if not 0 < length <= 16000:
                raise ValueError("Invalid request size")
            body = json.loads(self.rfile.read(length))
            self.json(run(body))
        except (ValueError, TypeError, KeyError) as e:
            self.json({"error": str(e)}, 400)
        finally:
            COMPUTE.release()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--port", type=int, default=8841)
    args = p.parse_args()
    load()
    frozen_model()
    http = ThreadingHTTPServer(("127.0.0.1", args.port), Handler)
    print(f"Spatial tumor lab: http://127.0.0.1:{http.server_port}", flush=True)
    try:
        http.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        http.server_close()


if __name__ == "__main__":
    main()
