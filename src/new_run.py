"""./new-run <run name>  --  create an empty run folder with a settings file to fill in."""

import re
import sys
from datetime import date

from runs import ROOT, RUNS, list_runs


def main():
    if len(sys.argv) != 2:
        print("Usage: ./new-run <run name>     e.g. ./new-run grafted_tierA_neutral_wall\n")
        list_runs()
        sys.exit(1)

    name = sys.argv[1]
    if not re.fullmatch(r"[A-Za-z0-9._-]+", name):
        sys.exit("A run name can only use letters, digits, '.', '_' and '-'.")
    run_dir = RUNS / name
    if run_dir.exists():
        sys.exit(f"runs/{name} already exists. Pick another name.")

    run_dir.mkdir(parents=True)
    template = (ROOT / "run_template.yaml").read_text()
    template = template.replace("RUN_NAME", name).replace("CREATED", date.today().isoformat())
    (run_dir / "run.yaml").write_text(template)

    print(f"Created runs/{name}/\n")
    print("Next:")
    print(f"  1. copy your sequence CSV into runs/{name}/")
    print(f"  2. edit runs/{name}/run.yaml")
    print(f"  3. ./launch-run {name}")


if __name__ == "__main__":
    main()
