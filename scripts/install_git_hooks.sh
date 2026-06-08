#!/usr/bin/env bash
# Install the repo's git hooks (privacy pre-commit guard).
source "$(dirname "$0")/lib/common.sh"
git -C "$REPO_ROOT" config core.hooksPath .githooks
chmod +x "$REPO_ROOT/.githooks/"* 2>/dev/null || true
ok "core.hooksPath set to .githooks — pre-commit privacy guard is active."
