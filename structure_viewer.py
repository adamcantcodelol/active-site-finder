"""Interactive 3D structure comparison helpers for MBRC Active Site Finder."""
from __future__ import annotations

import html
import json
import re

_PDB_RE = re.compile(r"^[0-9A-Za-z]{4}$")
_CHAIN_RE = re.compile(r"^[A-Za-z0-9_]+$")
_RES_RE = re.compile(r"^[+-]?\d+[A-Za-z]?$")


def normalize_pdb_id(pdb_id: str) -> str:
    raw = str(pdb_id).strip()
    if not raw:
        raise ValueError("PDB ID is empty")
    candidate = raw.split()[0][:4]
    if not _PDB_RE.fullmatch(candidate):
        raise ValueError("PDB ID must be exactly four letters/numbers")
    return candidate.lower()


def normalize_residues(residues) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for value in residues or []:
        text = str(value).strip()
        if _RES_RE.fullmatch(text) and text not in seen:
            seen.add(text)
            out.append(text)
    return out


def viewer_html(
    pdb_id: str,
    chain: str | None = None,
    residues=None,
    height: int = 560,
    homolog_pdb_id: str | None = None,
    homolog_chain: str | None = None,
    homolog_residues=None,
) -> str:
    """Create a self-contained 3Dmol viewer with four comparison modes."""
    query = normalize_pdb_id(pdb_id)
    qchain = str(chain or "").strip()
    if qchain and not _CHAIN_RE.fullmatch(qchain):
        raise ValueError("Invalid query chain identifier")
    qres = normalize_residues(residues)

    homolog = normalize_pdb_id(homolog_pdb_id) if homolog_pdb_id else ""
    hchain = str(homolog_chain or "").strip()
    if hchain and not _CHAIN_RE.fullmatch(hchain):
        raise ValueError("Invalid homolog chain identifier")
    hres = normalize_residues(homolog_residues)

    safe_height = max(420, min(int(height), 900))
    payload = {"query": query.upper(), "qchain": qchain, "qres": qres,
               "homolog": homolog.upper(), "hchain": hchain, "hres": hres}
    data = json.dumps(payload).replace("</", "<\\/")

    return f'''<!doctype html>
<html><head><meta charset="utf-8">
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>html,body{{margin:0;padding:0;background:#fff;font-family:Arial,Helvetica,sans-serif}}#viewer{{width:100%;height:{safe_height-52}px;position:relative}}#status{{height:52px;display:flex;align-items:center;padding:0 12px;box-sizing:border-box;color:#64748b;font-size:13px;border-top:1px solid #e2e8f0}}</style>
</head><body><div id="viewer"></div><div id="status">Loading structures…</div>
<script>
const DATA={data};
let viewer=null, qModel=null, hModel=null;
const ORANGE='orange', CYAN='cyan';
function atoms(model, chain, resis) {{
  if(!model || !chain || !resis || !resis.length) return [];
  const wanted=new Set(resis.map(String));
  return model.selectedAtoms({{chain:chain}}).filter(a=>wanted.has(String(a.resi)));
}}
function sub(a,b) {{ return {{x:a.x-b.x,y:a.y-b.y,z:a.z-b.z}}; }}
function add(a,b) {{ return {{x:a.x+b.x,y:a.y+b.y,z:a.z+b.z}}; }}
function mul(a,s) {{ return {{x:a.x*s,y:a.y*s,z:a.z*s}}; }}
function dot(a,b) {{ return a.x*b.x+a.y*b.y+a.z*b.z; }}
function cross(a,b) {{ return {{x:a.y*b.z-a.z*b.y,y:a.z*b.x-a.x*b.z,z:a.x*b.y-a.y*b.x}}; }}
function norm(a) {{ const n=Math.sqrt(dot(a,a))||1; return mul(a,1/n); }}
function basis(points) {{
  const o=points[0];
  const e1=norm(sub(points[1],o));
  const raw=sub(points[2],o);
  const e2=norm(sub(raw,mul(e1,dot(raw,e1))));
  const e3=norm(cross(e1,e2));
  return {{o:o,e1:e1,e2:e2,e3:e3}};
}}
function toBasis(p,b) {{ const v=sub(p,b.o); return {{x:dot(v,b.e1),y:dot(v,b.e2),z:dot(v,b.e3)}}; }}
function fromBasis(v,b) {{ return add(b.o,add(mul(b.e1,v.x),add(mul(b.e2,v.y),mul(b.e3,v.z)))); }}
function alphaFor(model,chain,res) {{
  return model.selectedAtoms({{chain:chain,resi:String(res)}}).find(a=>a.atom==='CA'||a.name==='CA') || null;
}}
function overlayBySite() {{
  if(!qModel||!hModel||DATA.qres.length<3||DATA.hres.length<3) return false;
  const qc=[],hc=[];
  for(let i=0;i<3;i++) {{
    const q=alphaFor(qModel,DATA.qchain,DATA.qres[i]);
    const h=alphaFor(hModel,DATA.hchain,DATA.hres[i]);
    if(!q||!h) return false;
    qc.push({{x:q.x,y:q.y,z:q.z}}); hc.push({{x:h.x,y:h.y,z:h.z}});
  }}
  const qb=basis(qc), hb=basis(hc);
  for(const a of hModel.atoms) {{
    const local=toBasis({{x:a.x,y:a.y,z:a.z}},hb);
    const p=fromBasis(local,qb); a.x=p.x; a.y=p.y; a.z=p.z;
  }}
  return true;
}}
function clearStyles() {{
  qModel.setStyle({{}},{{cartoon:{{color:'#cbd5e1'}}}});
  hModel.setStyle({{}},{{cartoon:{{color:'#93c5fd'}}}});
}}
function showMode(mode) {{
  if(!viewer||!qModel||!hModel) return;
  clearStyles();
  const qSel={{chain:DATA.qchain,resi:DATA.qres}};
  const hSel={{chain:DATA.hchain,resi:DATA.hres}};
  if(mode===2 || mode===4) {{ qModel.setStyle({{}},{{}}); hModel.setStyle({{}},{{}}); }}
  qModel.setStyle(qSel,{{stick:{{radius:0.24,colorscheme:'orangeCarbon'}},sphere:{{radius:0.38,color:ORANGE}}}});
  hModel.setStyle(hSel,{{stick:{{radius:0.24,colorscheme:'cyanCarbon'}},sphere:{{radius:0.38,color:CYAN}}}});
  if(mode===3 || mode===4) {{ overlayBySite(); }}
  viewer.zoomTo(); viewer.render();
  const labels={{1:'Mode 1 — full proteins + highlighted 3-residue sites',2:'Mode 2 — active-site residues only',3:'Mode 3 — overlaid proteins + highlighted sites',4:'Mode 4 — overlaid active-site residues only'}};
  document.getElementById('status').textContent=labels[mode];
}}
async function loadPdb(id) {{
  const r=await fetch('https://files.rcsb.org/download/'+encodeURIComponent(id)+'.pdb');
  if(!r.ok) throw new Error('Could not download '+id+' from RCSB');
  return await r.text();
}}
async function init() {{
 try {{
  viewer=$3Dmol.createViewer(document.getElementById('viewer'),{{backgroundColor:'white'}});
  const qp=await loadPdb(DATA.query); qModel=viewer.addModel(qp,'pdb');
  if(!DATA.homolog) {{ document.getElementById('status').textContent='Select a SPRITE homolog to compare structures.'; return; }}
  const hp=await loadPdb(DATA.homolog); hModel=viewer.addModel(hp,'pdb');
  clearStyles();
  qModel.setStyle({{chain:DATA.qchain,resi:DATA.qres}},{{stick:{{radius:0.24,colorscheme:'orangeCarbon'}},sphere:{{radius:0.38,color:ORANGE}}}});
  hModel.setStyle({{chain:DATA.hchain,resi:DATA.hres}},{{stick:{{radius:0.24,colorscheme:'cyanCarbon'}},sphere:{{radius:0.38,color:CYAN}}}});
  viewer.zoomTo(); viewer.render();
  window.applyMode=showMode;
  document.getElementById('status').textContent='Mode 1 — full proteins + highlighted 3-residue sites';
 }} catch(e) {{ document.getElementById('status').textContent='3D viewer error: '+e.message; }}
}}
init();
</script></body></html>'''
