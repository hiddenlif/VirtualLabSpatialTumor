"""Maintenance must not silently admit unrecorded scientific source changes."""

import json
import unittest
from unittest.mock import patch

from scripts import run_study


class SourceBindingTests(unittest.TestCase):
    def setUp(self):
        record = run_study.ROOT / "reports/study-v1/RUN.json"
        self.baseline = json.loads(record.read_text())["scientific_sources"]

    def test_recorded_cleanup_revision_is_accepted(self):
        self.assertEqual(
            run_study.verify_source_binding(self.baseline), "ruff-cleanup-v1"
        )

    def test_unrecorded_source_edit_is_rejected(self):
        changed = run_study.scientific_sources()
        changed["spatial_lab/analysis.py"] = "0" * 64
        with (
            patch.object(run_study, "scientific_sources", return_value=changed),
            self.assertRaises(ValueError),
        ):
            run_study.verify_source_binding(self.baseline)

    def test_unrelated_original_study_is_rejected(self):
        self.baseline["spatial_lab/analysis.py"] = "0" * 64
        with self.assertRaises(ValueError):
            run_study.verify_source_binding(self.baseline)
