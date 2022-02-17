<!--
# qbic-pipelines/rnadeseq pull request

Many thanks for contributing to qbic-pipelines/rnadeseq!

Please fill in the appropriate checklist below (delete whatever is not relevant).
These are the most common things requested on pull requests (PRs).

Remember that PRs should be made against the dev branch, unless you're preparing a pipeline release.

Learn more about contributing: [CONTRIBUTING.md](https://github.com/qbic-pipelines/rnadeseq/tree/master/.github/CONTRIBUTING.md)
-->
<!-- markdownlint-disable ul-indent -->

## PR checklist
//TODO: Change nf-core? In "If necessary, also...", this does not make much sense, maybe just delete?

- [ ] This comment contains a description of changes (with reason).
- [ ] If you've fixed a bug or added code that should be tested, add tests!
    - [ ] If you've added a new tool - have you followed the pipeline conventions in the [contribution docs](https://github.com/qbic-pipelines/rnadeseq/tree/master/.github/CONTRIBUTING.md)
    - [ ] If necessary, also make a PR on the qbic-pipelines/rnadeseq _branch_ on the [nf-core/test-datasets](https://github.com/nf-core/test-datasets) repository.
- [ ] Make sure your code lints (`nf-core lint`).
- [ ] Ensure the test suite passes (`nextflow run . -profile test,docker`).
- [ ] Usage Documentation in `docs/usage.md` is updated.
- [ ] Output Documentation in `docs/output.md` is updated.
- [ ] `CHANGELOG.md` is updated.
- [ ] `README.md` is updated (including new tool citations and authors/contributors).
