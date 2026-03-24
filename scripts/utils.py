import yaml
import os
import csv


def load_config(path="config.yaml"):
    with open(path) as f:
        return yaml.safe_load(f)


def file_exists_skip(path, label=""):
    if os.path.exists(path):
        print(f"[skip] {label or path}")
        return True
    return False


def append_metric(csv_path: str, row: dict):
    os.makedirs(os.path.dirname(csv_path) or ".", exist_ok=True)
    fieldnames = list(row.keys())

    write_header = not os.path.exists(csv_path) or os.stat(csv_path).st_size == 0

    existing_fieldnames = []
    if os.path.exists(csv_path) and os.stat(csv_path).st_size > 0:
        with open(csv_path, "r", newline="") as f:
            reader = csv.DictReader(f)
            existing_fieldnames = reader.fieldnames or []

    all_fieldnames = list(dict.fromkeys(existing_fieldnames + fieldnames))

    if set(all_fieldnames) != set(existing_fieldnames) and existing_fieldnames:
        rows = []
        with open(csv_path, "r", newline="") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        rows.append(row)
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=all_fieldnames, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
    else:
        with open(csv_path, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=all_fieldnames, extrasaction="ignore")
            if write_header:
                writer.writeheader()
            writer.writerow(row)
