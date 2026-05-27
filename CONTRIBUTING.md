# Contributing to cadcore

Thank you for your interest in cadcore!

## What we welcome

- Bug reports and reproducible test cases
- Documentation improvements and typo fixes
- Performance improvements with benchmarks
- New geometric primitives (curves, surfaces)
- Additional export formats (3MF, IGES, OBJ)
- Test coverage improvements

## What requires discussion first

- New public API surface (open an issue first)
- Changes to the B-Rep data model
- Commercial use cases — contact dmytroyatskovskiy@gmail.com

## How to contribute

1. Fork the repository
2. Create a branch: `git checkout -b fix/my-fix`
3. Make your changes
4. Run tests: `cargo test --workspace`
5. Run clippy: `cargo clippy --workspace -- -D warnings`
6. Open a pull request

## Code style

- `cargo fmt` before committing
- All public items must have doc comments (`///`)
- New features must include tests
- No unsafe code without a clear justification

## License

By contributing, you agree that your contributions will be licensed under the
[PolyForm Noncommercial License 1.0.0](LICENSE).
