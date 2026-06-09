# genome-toolkit MCP server

Serve your processed genome as a private knowledge source any AI agent can
query — over [Model Context Protocol](https://modelcontextprotocol.io).

```
run_all.sh ─► analysis/*.tsv/vcf ─► build_insights.py ─► genome_insights.json ─► server.py ─► your agent
                                        (all gitignored, local)              (stdio, no network)
```

> **Privacy.** The server speaks **stdio** — your agent launches it as a local
> subprocess; it never opens a port or sends data anywhere. The insights file
> (`analysis/genome_insights.json`) is gitignored. Only this code is shared.

## Setup (one command)

```bash
scripts/50_serve_mcp.sh        # venv + install + build insights + print config
```

That creates `mcp/.venv`, installs the MCP SDK, builds
`analysis/genome_insights.json`, and prints the exact config to paste into your
agent. Re-run with `--rebuild` after re-processing your data.

## Connect your agent

**Claude Code (CLI):**
```bash
claude mcp add genome -- /ABS/PATH/mcp/.venv/bin/python /ABS/PATH/mcp/server.py
```

**Claude Desktop** — add to `claude_desktop_config.json`:
```json
{
  "mcpServers": {
    "genome": {
      "command": "/ABS/PATH/mcp/.venv/bin/python",
      "args": ["/ABS/PATH/mcp/server.py"]
    }
  }
}
```

**Cursor** — same shape in `~/.cursor/mcp.json`. Any MCP client works; it just
needs the `command` + `args` above.

**Test before connecting:**
```bash
scripts/50_serve_mcp.sh --inspect      # opens the MCP Inspector UI (needs npx)
```

## What the agent gets

**Tools**
| Tool | Does |
|---|---|
| `get_summary()` | sample, reference, QC, inferred sex, finding counts |
| `get_pharmacogenomics(actionable_only?)` | predicted drug-metabolizer phenotypes (CPIC) |
| `check_drug(drug)` | which of your pharmacogenes affect a given drug + guidance |
| `get_traits(category?)` | genotypes at well-known trait variants |
| `search_clinvar(gene?, tier?, min_stars?)` | ClinVar findings you carry, by tier/gene/confidence |
| `lookup_variant(query)` | **live** lookup of any rsID / position / region in your VCF |
| `semantic_search(query, …)` | **meaning-based** search over your findings (needs the local index) |

**Resources:** `genome://summary`, `genome://pharmacogenomics`, `genome://traits`

### Semantic search (optional, fully local)

Exact tools (`check_drug`, `search_clinvar`) match keywords. For open-ended
questions where the words won't appear in your data — *"anything about
connective tissue?"*, *"heart rhythm"*, *"serotonin"* — build the local
semantic index:

```bash
scripts/50_serve_mcp.sh --embeddings    # installs fastembed, embeds findings ON-DEVICE
```

Embeddings run on your CPU (fastembed/ONNX, no torch). Only a small **public**
model (~130 MB) is downloaded once; **your genome text is never sent anywhere**.
The index (`analysis/genome_embeddings.*`) is gitignored. The `semantic_search`
tool then ranks your findings by meaning. Rebuild it after re-processing data.

> **Design note — where the "smarts" live.** The agent handles *world* semantics
> (turning "is this supplement risky?" into "check CYP3A4 / serotonin genes" —
> knowledge your genome doesn't contain). This server handles *data* semantics
> (matching that concept to your rows even when the words differ). Keep
> world-knowledge in the agent; keep your-genome retrieval here.

## Example questions your agent can now answer

- "I might be prescribed clopidogrel — any genetic concern?" → `check_drug("clopidogrel")`
- "What's my APOE status?" → `lookup_variant("rs429358")` + `lookup_variant("rs7412")`
- "Summarize my pharmacogenomic risks." → `get_pharmacogenomics(actionable_only=True)`
- "Do I carry anything pathogenic in DYSF?" → `search_clinvar(gene="DYSF")`

## Remote hosting (advanced, discouraged)

This is genomic data — keep it local. If you genuinely need a remote server
(e.g. a shared agent), do **not** expose it unauthenticated. The MCP SDK
supports a streamable-HTTP transport; put it behind a reverse proxy with auth
and TLS, bind to localhost + an SSH tunnel, or use your client's remote-MCP
gateway. Treat the endpoint like a password to your medical record.

## ⚠️ Disclaimer

Educational only — not medical advice. Predictions come from consumer-grade WGS.
Confirm anything actionable with a clinician or pharmacist (pharmacogenomics:
[CPIC](https://cpicpgx.org)).
