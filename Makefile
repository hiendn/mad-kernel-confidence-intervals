.PHONY: check verify
verify:
	Rscript scripts/verify_workspace.R
check:
	Rscript scripts/check.R
