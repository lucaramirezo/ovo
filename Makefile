# Ovo: one entry point for the CAD build, media, report and site.
#
#   make setup    Python env (.venv) from uv.lock, docs group included
#   make build    CAD: part modules, assembly checks, STEP / STL / GLB, renders (scripts/build.py)
#   make fresh    the same from scratch: drop the BREP caches, prune stale outputs (as CI does)
#   make media    README / site media: hero stills, fold GIF and MP4 (scripts/media.py)
#   make report   technical report PDF, docs/report/ovo-report.pdf (scripts/report.py)
#   make site     project web page assets under site/ (scripts/site.py)
#   make all      build, media, report, site in that order (each later step reads the earlier ones)
#   make links    check that the repository, release, Pages and og:image links answer (before sharing)
#
# Tools: uv (https://docs.astral.sh/uv/) and, for the viewer check and every render,
# agent-browser with its headless Chrome: npm install -g agent-browser && agent-browser install
# Heavy steps take out/work/heavy.lock and scripts/cadrun.sh caps OCC jobs at 2, so parallel
# make jobs would only queue: the file runs serially.

PY  ?= .venv/bin/python
UV  ?= uv
BUILD_FLAGS ?=
PORT ?= 8000

.DEFAULT_GOAL := help
.NOTPARALLEL:
.PHONY: help setup build fresh summary media report site all links serve-site viewer check-tools

help:
	@sed -n '3,10p' $(firstword $(MAKEFILE_LIST)) | sed 's/^# \{0,1\}//'

setup:
	$(UV) sync --locked --all-groups
	@command -v agent-browser > /dev/null 2>&1 || \
	  echo "note: the build and the media need agent-browser: npm install -g agent-browser && agent-browser install"

check-tools:
	@test -x $(PY) || { echo "no $(PY): run 'make setup' first"; exit 1; }
	@command -v agent-browser > /dev/null 2>&1 || \
	  { echo "agent-browser not found: npm install -g agent-browser && agent-browser install"; exit 1; }
	@command -v flock > /dev/null 2>&1 || \
	  { echo "flock not found (util-linux; macOS: brew install flock): scripts/cadrun.sh needs it"; exit 1; }

build: check-tools
	$(PY) scripts/build.py $(BUILD_FLAGS)

fresh: check-tools
	$(PY) scripts/build.py --fresh $(BUILD_FLAGS)

summary:
	$(PY) scripts/build.py --summary

# media, report and site read the build outputs; they trigger a build only when those are missing
out/web/ovo.glb out/cad/checks/assembly.json:
	$(MAKE) build

media: check-tools | out/web/ovo.glb
	$(PY) scripts/media.py

report: | out/cad/checks/assembly.json
	$(PY) scripts/report.py

# site/index.html is hand-written; scripts/site.py refreshes its assets (compressed GLB via npx
# @gltf-transform/cli, manifest, media and report copies, generated data blocks) and checks its numbers
site: | out/web/ovo.glb
	$(PY) scripts/site.py

all: build media report site

links:
	$(PY) scripts/site.py --links

# Local preview servers (ES modules and the GLB need http, not file://)
serve-site:
	$(PY) -m http.server $(PORT) -d site

viewer:
	@echo "viewer: http://127.0.0.1:$(PORT)/viewer/index.html?ui=1  (URL parameters: the comment at the top of its script)"
	$(PY) -m http.server $(PORT)
