# Minimal Makefile for the tutorial site.
# Point these at your sphinx-docs env (e.g. SPHINXBUILD=/path/to/conda_envs/sphinx-docs/bin/sphinx-build),
# or activate the env so sphinx-build / sphinx-autobuild are on PATH.
SPHINXBUILD ?= sphinx-build
AUTOBUILD   ?= sphinx-autobuild
SOURCEDIR    = source
BUILDDIR     = build

html:
	$(SPHINXBUILD) -b html "$(SOURCEDIR)" "$(BUILDDIR)/html"

strict:                       ## fail the build on broken xrefs / missing images
	$(SPHINXBUILD) -W -b html "$(SOURCEDIR)" "$(BUILDDIR)/html"

serve:                        ## live-reloading local preview
	$(AUTOBUILD) "$(SOURCEDIR)" "$(BUILDDIR)/html"

clean:
	rm -rf "$(BUILDDIR)"

.PHONY: html strict serve clean
