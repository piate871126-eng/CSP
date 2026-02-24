from __future__ import annotations

import base64
import json
from pathlib import Path
from uuid import uuid4

from flask import Flask, render_template, request

from .pipeline import run_pipeline_from_smiles
from .smiles import SmilesConversionError


def create_app() -> Flask:
    app = Flask(__name__, template_folder="templates")

    @app.route("/", methods=["GET", "POST"])
    def index():
        ctx: dict = {"smiles": "", "error": None, "report": None, "xrd_image_b64": None}
        if request.method == "POST":
            smiles = request.form.get("smiles", "").strip()
            ctx["smiles"] = smiles
            if not smiles:
                ctx["error"] = "SMILES를 입력해 주세요."
                return render_template("index.html", **ctx)

            out_dir = Path("outputs") / f"web_{uuid4().hex[:8]}"
            try:
                report_path = run_pipeline_from_smiles(smiles, out_dir)
                report = json.loads(report_path.read_text(encoding="utf-8"))
                img_path = Path(report["xrd"]["plot"])
                xrd_image_b64 = base64.b64encode(img_path.read_bytes()).decode("ascii")
                ctx["report"] = report
                ctx["xrd_image_b64"] = xrd_image_b64
            except SmilesConversionError as exc:
                ctx["error"] = str(exc)
            except Exception as exc:  # pragma: no cover
                ctx["error"] = f"처리 중 오류가 발생했습니다: {exc}"

        return render_template("index.html", **ctx)

    return app


def main() -> None:
    app = create_app()
    app.run(host="0.0.0.0", port=8000, debug=False)


if __name__ == "__main__":
    main()
