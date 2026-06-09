#!/usr/bin/env bash
# ============================================================================
# 50_serve_mcp.sh — set up and serve your genome insights over MCP.
#
# One command: creates a venv, installs the MCP SDK, builds the structured
# insights file, and prints the exact config to paste into your agent
# (Claude Desktop / Claude Code / Cursor). Everything stays local — the server
# speaks stdio and never opens a network port.
#
#   scripts/50_serve_mcp.sh            # setup + build + print connection config
#   scripts/50_serve_mcp.sh --rebuild  # just rebuild genome_insights.json
#   scripts/50_serve_mcp.sh --run      # run the server in the foreground (stdio)
#   scripts/50_serve_mcp.sh --inspect  # launch the MCP Inspector UI to test it
# ============================================================================
source "$(dirname "$0")/lib/common.sh"
require_tool bcftools

MCP_DIR="$REPO_ROOT/mcp"
VENV="$MCP_DIR/.venv"
PYBIN="$VENV/bin/python"
INSIGHTS="$OUT_DIR/genome_insights.json"

ensure_venv() {
  if [[ ! -x "$PYBIN" ]]; then
    log "Creating venv at $VENV…"
    python3 -m venv "$VENV"
    "$PYBIN" -m pip install -q --upgrade pip
    "$PYBIN" -m pip install -q -r "$MCP_DIR/requirements.txt"
    ok "Installed MCP SDK"
  fi
}

build_insights() {
  [[ -f "$OUT_DIR/variants.vcf.gz" ]] || die "No variants.vcf.gz — run scripts/run_all.sh first."
  log "Building structured insights…"
  OUT_DIR="$OUT_DIR" "$PYBIN" "$MCP_DIR/build_insights.py"
}

print_config() {
  hr
  log "MCP server is ready. Connect it to your agent:"
  cat >&2 <<EOF

  Command + args (any MCP client):
      command : $PYBIN
      args    : ["$MCP_DIR/server.py"]

  ── Claude Code (CLI) ──────────────────────────────────────────────
      claude mcp add genome -- "$PYBIN" "$MCP_DIR/server.py"

  ── Claude Desktop ── add to claude_desktop_config.json:
  {
    "mcpServers": {
      "genome": {
        "command": "$PYBIN",
        "args": ["$MCP_DIR/server.py"]
      }
    }
  }

  ── Cursor ── add to ~/.cursor/mcp.json (same shape as above).

  Insights file served: $INSIGHTS
  Test it first with:   scripts/50_serve_mcp.sh --inspect
EOF
  hr
}

case "${1:-}" in
  --rebuild) ensure_venv; build_insights ;;
  --run)     ensure_venv; [[ -f "$INSIGHTS" ]] || build_insights; exec "$PYBIN" "$MCP_DIR/server.py" ;;
  --inspect) ensure_venv; [[ -f "$INSIGHTS" ]] || build_insights
             log "Launching MCP Inspector (needs npx)…"
             exec npx -y @modelcontextprotocol/inspector "$PYBIN" "$MCP_DIR/server.py" ;;
  *)         ensure_venv; build_insights; print_config ;;
esac
