#!/usr/bin/env bash
set -euo pipefail

repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root"

for required_command in uv cmake ninja c++ sha256sum; do
    if ! command -v "$required_command" >/dev/null 2>&1; then
        echo "missing required command: $required_command" >&2
        exit 1
    fi
done

git_common_dir="$(git rev-parse --path-format=absolute --git-common-dir)"
canonical_repo_root="$(dirname "$git_common_dir")"
project_root="$(dirname "$canonical_repo_root")"
default_eos_wheel="${project_root}/artifacts/objects/sha256/66b7ea8fb29e0a268b555cbdf401c3502517c088669a4157e8f64ab985b59ce9/epcsaft-0.2.0.dev0-cp313-cp313-linux_x86_64.whl"
eos_wheel="${EPCSAFT_EOS_WHEEL:-$default_eos_wheel}"
expected_eos_sha256="66b7ea8fb29e0a268b555cbdf401c3502517c088669a4157e8f64ab985b59ce9"

if [[ ! -f "$eos_wheel" ]]; then
    echo "missing required immutable EOS wheel: $eos_wheel" >&2
    echo "set EPCSAFT_EOS_WHEEL to the exact retained wheel if the artifact root moved" >&2
    exit 1
fi

actual_eos_sha256="$(sha256sum "$eos_wheel" | cut -d ' ' -f 1)"
if [[ "$actual_eos_sha256" != "$expected_eos_sha256" ]]; then
    echo "EOS wheel SHA-256 mismatch: expected $expected_eos_sha256, got $actual_eos_sha256" >&2
    exit 1
fi

if [[ ! -x .venv/bin/python ]]; then
    uv venv --python 3.13 .venv
fi

uv pip install --python .venv/bin/python --reinstall-package epcsaft \
    "$eos_wheel" \
    pytest \
    scikit-build-core

eos_include_dir="$(
    .venv/bin/python - <<'PY'
from pathlib import Path
import epcsaft

print(Path(epcsaft.__file__).resolve().parent / "include")
PY
)"
export EPCSAFT_INCLUDE_DIR="$eos_include_dir"

uv pip install --python .venv/bin/python --no-build-isolation --editable .

.venv/bin/python - "$eos_wheel" <<'PY'
from importlib import metadata
import json
from pathlib import Path
import sys
from urllib.parse import unquote, urlparse

expected_wheel = Path(sys.argv[1]).resolve()
eos = metadata.distribution("epcsaft")
direct_url = json.loads(eos.read_text("direct_url.json") or "{}")
if eos.version != "0.2.0.dev0":
    raise SystemExit(f"epcsaft version mismatch: {eos.version}")
if direct_url.get("dir_info", {}).get("editable", False):
    raise SystemExit("editable epcsaft EOS install is forbidden")
installed_url = direct_url.get("url", "")
installed_wheel = Path(unquote(urlparse(installed_url).path)).resolve()
if installed_wheel != expected_wheel:
    raise SystemExit(
        f"installed EOS wheel mismatch: {installed_wheel} != {expected_wheel}"
    )
if "archive_info" not in direct_url:
    raise SystemExit("epcsaft EOS must be installed from an immutable wheel")
PY

cmake -S . -B build -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DPython_EXECUTABLE="$repo_root/.venv/bin/python"
cmake --build build -j2
