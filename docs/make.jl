# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

using Documenter
using BoltzTraP

makedocs(
    sitename = "BoltzTraP.jl",
    authors = "Hiroharu Sugawara",
    modules = [BoltzTraP],
    remotes = nothing,  # Disable source links when not in Git repo
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://hsugawa8651.github.io/BoltzTraP.jl",
        edit_link = "main",
    ),
    pages = [
        "Home" => "index.md",
        "Guide" => [
            "CLI Workflow" => "cli.md",
            "API Workflow" => "workflow.md",
            "DFTK.jl Integration" => "dftk.md",
            "Interpolation" => "interpolate.md",
            "Integration" => "integrate.md",
            "Plotting" => "plotting.md",
            "Input Formats" => "input_formats.md",
            "Output Formats" => "output_formats.md",
            "Conventions" => "conventions.md",
            "Functional Tests" => "ftest.md",
        ],
        "Validation" => "validation.md",
        "Benchmarks" => "benchmarks.md",
        "Reference Tests" => "reftest.md",
    ],
    warnonly = [:missing_docs, :cross_references],
)

deploydocs(
    repo = "github.com/hsugawa8651/BoltzTraP.jl.git",
    devbranch = "main",
    push_preview = true,
)
