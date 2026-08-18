# ACG-BRNS Documentation

Diese Dokumentation wird mit MkDocs und mkdocstrings erzeugt.

## Inhalt

- API-Dokumentation fuer die Kernmodule:
  - `acg_brns.acg`
  - `acg_brns.gaussian_elimination`
  - `acg_brns.acg_orchestrator`

## Lokaler Build

```bash
pip install -e .
pip install mkdocs-material "mkdocstrings[python]" pymdown-extensions
mkdocs serve
```

Danach ist die Doku lokal unter `http://127.0.0.1:8000` verfuegbar.
