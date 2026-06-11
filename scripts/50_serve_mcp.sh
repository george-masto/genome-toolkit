#!/usr/bin/env bash
# ============================================================================
# 50_serve_mcp.sh - set up and serve your genome insights over MCP.
#
# One command: creates a venv, installs the MCP SDK, builds the structured
# insights file, and prints the exact config to paste into your agent
# (Claude Desktop / Claude Code / Cursor). Everything stays local - the server
# speaks stdio and never opens a network port.
#
#   scripts/50_serve_mcp.sh               # setup + build + print stdio config
#   scripts/50_serve_mcp.sh --rebuild     # just rebuild genome_insights.json
#   scripts/50_serve_mcp.sh --embeddings  # build the LOCAL semantic-search index
#   scripts/50_serve_mcp.sh --run         # run in the foreground (stdio)
#   scripts/50_serve_mcp.sh --inspect     # launch the MCP Inspector UI to test it
#   scripts/50_serve_mcp.sh --http        # serve over HTTP bound to your Tailscale IP
#   scripts/50_serve_mcp.sh --install-service    # always-on launchd service (HTTP/Tailscale)
#   scripts/50_serve_mcp.sh --uninstall-service  # remove the launchd service
# ============================================================================
source "$(dirname "$0")/lib/common.sh"
require_tool bcftools

MCP_DIR="$REPO_ROOT/mcp"
VENV="$MCP_DIR/.venv"
PYBIN="$VENV/bin/python"
INSIGHTS="$OUT_DIR/genome_insights.json"

ensure_venv() {
  if [[ ! -x "$PYBIN" ]]; then
    log "Creating venv at ${VENV}..."
    "${PYTHON:-python3}" -m venv "$VENV"
    "$PYBIN" -m pip install -q --upgrade pip
    "$PYBIN" -m pip install -q -r "$MCP_DIR/requirements.txt"
    ok "Installed MCP SDK"
  fi
}

build_insights() {
  [[ -f "$OUT_DIR/variants.vcf.gz" ]] || die "No variants.vcf.gz - run scripts/run_all.sh first."
  log "Building structured insights..."
  OUT_DIR="$OUT_DIR" "$PYBIN" "$MCP_DIR/build_insights.py"
}

build_embeddings() {
  [[ -f "$INSIGHTS" ]] || build_insights
  log "Installing local embedding deps (one-time)..."
  "$PYBIN" -m pip install -q -r "$MCP_DIR/requirements-embeddings.txt"
  log "Building LOCAL semantic-search index (model downloads once, then offline)..."
  OUT_DIR="$OUT_DIR" "$PYBIN" "$MCP_DIR/build_embeddings.py"
  ok "Semantic search ready - the server's semantic_search tool is now live."
}

TOKEN_FILE="$MCP_DIR/.mcp_token"
PLIST="$HOME/Library/LaunchAgents/com.genome-toolkit.mcp.plist"

tailscale_ip() { tailscale ip -4 2>/dev/null | head -1; }

ensure_token() {
  if [[ -n "${MCP_AUTH_TOKEN:-}" ]]; then printf '%s' "$MCP_AUTH_TOKEN"; return; fi
  [[ -f "$TOKEN_FILE" ]] || { openssl rand -hex 24 > "$TOKEN_FILE"; chmod 600 "$TOKEN_FILE"; }
  cat "$TOKEN_FILE"
}

http_url() {  # echoes host port token via globals
  HOST="${MCP_HOST:-$(tailscale_ip)}"; HOST="${HOST:-127.0.0.1}"
  PORT="${MCP_PORT:-8787}"
  TOKEN="$(ensure_token)"
}

print_http_config() {
  http_url
  local ts_name; ts_name="$(tailscale status --json 2>/dev/null | sed -n 's/.*"DNSName":"\([^"]*\)".*/\1/p' | head -1 | sed 's/\.$//')"
  hr
  log "Genome MCP serving over HTTP on your tailnet:"
  cat >&2 <<EOF

  URL:    http://$HOST:$PORT/mcp
  ${ts_name:+MagicDNS URL: http://$ts_name:$PORT/mcp}
  Header: Authorization: Bearer $TOKEN

  ── Claude Code (CLI), from any tailnet machine ──────────────────
      claude mcp add --transport http genome http://$HOST:$PORT/mcp \\
        --header "Authorization: Bearer $TOKEN"

  ── Atlas / Claude Desktop (remote MCP / custom connector) ───────
      URL:     http://$HOST:$PORT/mcp
      Header:  Authorization = Bearer $TOKEN

  (Token stored at $TOKEN_FILE - gitignored. Reachable only inside your tailnet.)
EOF
  hr
}

serve_http() {
  ensure_venv
  [[ -f "$INSIGHTS" ]] || build_insights
  http_url
  print_http_config
  exec env MCP_TRANSPORT=http MCP_HOST="$HOST" MCP_PORT="$PORT" MCP_AUTH_TOKEN="$TOKEN" \
       GENOME_INSIGHTS="$INSIGHTS" "$PYBIN" "$MCP_DIR/server.py"
}

install_service() {
  ensure_venv
  [[ -f "$INSIGHTS" ]] || build_insights
  have tailscale || warn "tailscale CLI not found - the service will bind to 127.0.0.1 until it is."
  mkdir -p "$HOME/Library/LaunchAgents" "$REPO_ROOT/logs"
  local token; token="$(ensure_token)"
  cat > "$PLIST" <<EOF
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
  <key>Label</key><string>com.genome-toolkit.mcp</string>
  <key>ProgramArguments</key>
  <array>
    <string>/bin/bash</string>
    <string>$REPO_ROOT/scripts/50_serve_mcp.sh</string>
    <string>--http</string>
  </array>
  <key>EnvironmentVariables</key>
  <dict>
    <key>PATH</key><string>/opt/homebrew/bin:/usr/local/bin:/usr/bin:/bin</string>
    <key>MCP_AUTH_TOKEN</key><string>$token</string>
    <key>MCP_PORT</key><string>${MCP_PORT:-8787}</string>
  </dict>
  <key>RunAtLoad</key><true/>
  <key>KeepAlive</key><true/>
  <key>StandardOutPath</key><string>$REPO_ROOT/logs/mcp.out.log</string>
  <key>StandardErrorPath</key><string>$REPO_ROOT/logs/mcp.err.log</string>
</dict>
</plist>
EOF
  launchctl unload "$PLIST" 2>/dev/null || true
  launchctl load "$PLIST"
  ok "Installed always-on service: $PLIST"
  ok "Logs: $REPO_ROOT/logs/mcp.{out,err}.log"
  sleep 2; print_http_config
}

uninstall_service() {
  launchctl unload "$PLIST" 2>/dev/null || true
  rm -f "$PLIST"
  ok "Removed launchd service. (Token + data untouched.)"
}

# The server reads genome_insights.json once at startup, so a rebuild only
# reaches the LIVE MCP after a restart. If the always-on service is loaded,
# bounce it so a fresh build (e.g. a new CYP2D6 call) is served immediately.
restart_service_if_running() {
  if launchctl list 2>/dev/null | grep -q 'com\.genome-toolkit\.mcp'; then
    log "Restarting the always-on MCP service so it serves the rebuilt insights..."
    launchctl kickstart -k "gui/$(id -u)/com.genome-toolkit.mcp" 2>/dev/null \
      || { launchctl unload "$PLIST" 2>/dev/null; launchctl load "$PLIST" 2>/dev/null; }
    ok "MCP service restarted - live data is now the rebuilt insights."
  else
    warn "Insights rebuilt. No always-on service loaded; restart your MCP client/server to serve them."
  fi
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
  --rebuild)    ensure_venv; build_insights; restart_service_if_running ;;
  --embeddings) ensure_venv; build_embeddings ;;
  --http)       serve_http ;;
  --install-service)   install_service ;;
  --uninstall-service) uninstall_service ;;
  --run)     ensure_venv; [[ -f "$INSIGHTS" ]] || build_insights; exec "$PYBIN" "$MCP_DIR/server.py" ;;
  --inspect) ensure_venv; [[ -f "$INSIGHTS" ]] || build_insights
             log "Launching MCP Inspector (needs npx)..."
             exec npx -y @modelcontextprotocol/inspector "$PYBIN" "$MCP_DIR/server.py" ;;
  *)         ensure_venv; build_insights; print_config ;;
esac
