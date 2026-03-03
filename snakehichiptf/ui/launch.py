#!/usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess
from importlib import resources

def main():
    # robust way to locate packaged file
    app_path = resources.files("snakehichiptf.ui").joinpath("stl_app.py")
    cmd = [
        "streamlit", "run", str(app_path),
        "--browser.gatherUsageStats", "false",
    ]
    raise SystemExit(subprocess.call(cmd))

if __name__ == "__main__":
    main()