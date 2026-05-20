from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


def ensure_dir(path: str | Path) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def read_table(path: str | Path, sep: str | None = None, header: int | str | None = "infer") -> pd.DataFrame:
    path = Path(path)
    if sep is None:
        sep = "\t" if path.suffix.lower() in {".bed", ".txt", ".tsv"} else ","
    return pd.read_csv(path, sep=sep, header=header)


def maybe_named_columns(df: pd.DataFrame, names: Iterable[str]) -> pd.DataFrame:
    names = list(names)
    if all(isinstance(c, str) for c in df.columns) and set(names).issubset(df.columns):
        return df
    out = df.copy()
    out.columns = names[: len(out.columns)]
    return out


def write_json(path: str | Path, data: dict) -> None:
    with Path(path).open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2)


def read_json(path: str | Path) -> dict:
    with Path(path).open("r", encoding="utf-8") as handle:
        return json.load(handle)


def optional_pyplot():
    try:
        import matplotlib.pyplot as plt  # type: ignore

        return plt
    except Exception:
        return None


def interp1(x: np.ndarray, y: np.ndarray, xq: np.ndarray) -> np.ndarray:
    order = np.argsort(x)
    return np.interp(xq, np.asarray(x)[order], np.asarray(y)[order])


def save_data_struct(path_prefix: str | Path, **arrays: np.ndarray) -> None:
    path_prefix = Path(path_prefix)
    ensure_dir(path_prefix.parent)
    np.savez(path_prefix.with_suffix(".npz"), **arrays)
    write_json(path_prefix.with_suffix(".json"), {k: np.asarray(v).tolist() for k, v in arrays.items()})


def load_data_struct(path: str | Path) -> dict[str, np.ndarray]:
    path = Path(path)
    if path.suffix == ".npz":
        data = np.load(path)
        return {k: np.asarray(data[k], dtype=float) for k in data.files}
    if path.suffix == ".json":
        raw = read_json(path)
        return {k: np.asarray(v, dtype=float) for k, v in raw.items()}
    if path.suffix == ".mat":
        try:
            from scipy.io import loadmat  # type: ignore
        except Exception as exc:
            raise RuntimeError("Reading MATLAB .mat files requires scipy. Use the .npz/.json export instead.") from exc
        raw = loadmat(path, squeeze_me=True, struct_as_record=False)
        ds = raw.get("DataStruct")
        if ds is None:
            return {k: np.asarray(v, dtype=float).squeeze() for k, v in raw.items() if not k.startswith("__")}
        return {field: np.asarray(getattr(ds, field), dtype=float).squeeze() for field in ds._fieldnames}
    raise ValueError(f"Unsupported data-struct format: {path}")

