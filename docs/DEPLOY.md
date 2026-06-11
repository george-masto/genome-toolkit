# Deploying the genome MCP on your always-on server (Tailscale)

Run the genome MCP as an always-on HTTP service on your personal server, reachable
**only inside your tailnet** (your private network), so the AI assistants on your
other devices — Claude Desktop, Claude Code, and similar — can all use it.

```
 your laptop ───┐
 your phone  ───┤── Tailscale (WireGuard) ──► always-on server : http://<ts-ip>:8787/mcp
 Claude Desktop ┘                              genome-toolkit MCP (launchd, token-authed)
```

> **Privacy model.** The endpoint binds to the machine's **Tailscale IP**, so only
> tailnet devices can reach it; traffic is WireGuard-encrypted; a bearer token adds
> a second layer. Your genome data lives only on the server (gitignored). Nothing is
> public, nothing is committed.

---

## Step 1 — code on the server

```bash
ssh <you>@<server>            # or: tailscale ssh <server>
git clone https://github.com/george-masto/genome-toolkit ~/genome-toolkit
cd ~/genome-toolkit
scripts/00_install_tools.sh   # bcftools/samtools/htslib (Homebrew)
```

## Step 2 — get your data onto the server

The repo is code only — your genome is never in git, so transfer it directly over
the encrypted tailnet. Two options:

**Option A — self-contained server (recommended).** Send your raw download and
re-derive everything on the server (fully reproducible; needs ~120 GB working space):

```bash
# from THIS machine — point at your local data, then push it over Tailscale (encrypted)
LOCAL=/path/to/your/genome-data        # the folder holding your *.fq.gz / *.vcf.gz
rsync -avP "$LOCAL"/*.fq.gz "$LOCAL"/*.vcf.gz <you>@<server>:~/genome-toolkit/
# then ON THE SERVER, build everything:
scripts/run_all.sh                 # variants, QC, traits, drug-response
scripts/20_clinvar_annotate.sh     # ClinVar health-variant screen
```

**Option B — MCP-only (light, ~200 MB, no raw reads).** Copy just the processed
files the server needs:

```bash
LOCAL=/path/to/your/genome-toolkit     # your local checkout (results live in analysis/)
rsync -avP \
  "$LOCAL"/analysis/variants.vcf.gz \
  "$LOCAL"/analysis/variants.vcf.gz.tbi \
  "$LOCAL"/analysis/genome_insights.json \
  "$LOCAL"/analysis/genome_embeddings.npz \
  "$LOCAL"/analysis/genome_embeddings.jsonl \
  "$LOCAL"/analysis/20_clinvar_all_categorized.tsv \
  <you>@<server>:~/genome-toolkit/analysis/
```

(`rsync` runs over SSH; enable Remote Login on the server or use `tailscale ssh`.)

## Step 3 — build the index + semantic search

```bash
# ON THE SERVER
scripts/50_serve_mcp.sh --embeddings   # venv + insights (if needed) + local embeddings
```

For Option B this just installs deps and rebuilds the (already-copied) index — fast.

## Step 4 — run it always-on

```bash
scripts/50_serve_mcp.sh --install-service
```

This installs a **launchd** agent (`~/Library/LaunchAgents/com.genome-toolkit.mcp.plist`)
that, on login/boot and forever after (`KeepAlive`), runs the HTTP server bound to your
Tailscale IP with a token. It prints the **URL + token** to paste into your agents.
Logs: `~/genome-toolkit/logs/mcp.{out,err}.log`. (Foreground test first with
`scripts/50_serve_mcp.sh --http`.)

## Step 5 — connect your agents (from any tailnet device)

The setup command prints these with your real IP/token:

**Claude Code (CLI):**
```bash
claude mcp add --transport http genome http://<ts-ip>:8787/mcp \
  --header "Authorization: Bearer <token>"
```

**Claude Desktop / other MCP client** — add a remote MCP / custom connector:
- URL: `http://<ts-ip>:8787/mcp`
- Header: `Authorization: Bearer <token>`

MagicDNS works too: `http://<machine>.<tailnet>.ts.net:8787/mcp`.

---

## Maintenance

- **Re-process after new data/analysis:** rebuild and restart the service:
  ```bash
  scripts/50_serve_mcp.sh --rebuild        # genome_insights.json
  scripts/50_serve_mcp.sh --embeddings     # semantic index
  launchctl kickstart -k gui/$(id -u)/com.genome-toolkit.mcp   # restart service
  ```
- **Rotate the token:** `rm mcp/.mcp_token`, re-run `--install-service`, update agents.
- **Stop / remove:** `scripts/50_serve_mcp.sh --uninstall-service`.
- **Pin allowed hosts** (instead of disabling the Host check): set
  `MCP_ALLOWED_HOSTS="<ts-ip>:*,<machine>.<tailnet>.ts.net:*"` in the plist env.

## Hardening notes

- The network boundary (tailnet) is the primary control; the token is defense in depth.
- Keep the token out of git — it lives in `mcp/.mcp_token` (gitignored).
- Don't bind to `0.0.0.0` (that exposes it on the LAN too). `--http` binds to the
  Tailscale IP precisely to avoid that.
- This serves real health data — treat the URL+token like access to a medical record.
