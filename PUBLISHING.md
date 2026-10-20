# Publishing this portfolio repository

This directory is designed to become an independent GitHub repository.

From inside `DataDrivenFlowROM.jl`:

```bash
git init
git add .
git commit -m "Initial public release of DataDrivenFlowROM.jl"
git branch -M main
git remote add origin https://github.com/franabaajo/DataDrivenFlowROM.jl.git
git push -u origin main
```

Before publishing, run the test suite on a machine with Julia 1.10+ and check
that the CI badge resolves after the first GitHub Actions run. The thesis PDF is
about 7.5 MB; if repository size becomes a concern, attach it to a GitHub
release and replace the README link with the release URL.

