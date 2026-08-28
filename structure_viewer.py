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
    mode: int = 1,
) -> str:
    """Create a self-contained 3Dmol viewer with four SPRITE comparison modes."""
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

    if homolog and (len(qres) != 3 or len(hres) != 3):
        raise ValueError("A SPRITE comparison requires exactly three residues on each structure")
    if mode not in (1, 2, 3, 4):
        mode = 1

    safe_height = max(420, min(int(height), 900))
    payload = {"query": query.upper(), "qchain": qchain, "qres": qres,
               "homolog": homolog.upper(), "hchain": hchain, "hres": hres,
               "mode": mode}
    data = json.dumps(payload).replace("</", "<\\/")

    return f'''<!doctype html>
<html><head><meta charset="utf-8">
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>html,body{{margin:0;padding:0;background:#fff;font-family:Arial,Helvetica,sans-serif}}#viewer{{width:100%;height:{safe_height-52}px;position:relative}}#status{{height:52px;display:flex;align-items:center;padding:0 12px;box-sizing:border-box;color:#64748b;font-size:13px;border-top:1px solid #e2e8f0}}</style>
</head><body><div id="viewer"></div><div id="status">Loading structures…</div>
<script>
const DATA={data};
let viewer=null, qModel=null, hModel=null;
const BLUE='#3b82f6', YELLOW='#facc15', SITE_BLUE='#2563eb', SITE_YELLOW='#d97706';
function unique3(arr) {{
  const out=[]; const seen=new Set();
  for(const x of (arr||[])) {{ const s=String(x); if(!seen.has(s)) {{seen.add(s);out.push(s);}} }}
  return out.length===3 ? out : null;
}}
function ca(model, chain, res) {{
  const atoms=model.selectedAtoms({{chain:chain,resi:String(res)}});
  return atoms.find(a=>String(a.atom||'').toUpperCase()==='CA' || String(a.name||'').toUpperCase()==='CA') || null;
}}
function vec(a,b) {{return {{x:a.x-b.x,y:a.y-b.y,z:a.z-b.z}};}}
function add(a,b) {{return {{x:a.x+b.x,y:a.y+b.y,z:a.z+b.z}};}}
function scale(a,s) {{return {{x:a.x*s,y:a.y*s,z:a.z*s}};}}
function dot(a,b) {{return a.x*b.x+a.y*b.y+a.z*b.z;}}
function cross(a,b) {{return {{x:a.y*b.z-a.z*b.y,y:a.z*b.x-a.x*b.z,z:a.x*b.y-a.y*b.x}};}}
function unit(a) {{const n=Math.hypot(a.x,a.y,a.z); if(n<1e-8) throw new Error('SPRITE residues are geometrically degenerate'); return scale(a,1/n);}}
function frame(points) {{
  const origin=points[0];
  const e1=unit(vec(points[1],origin));
  const raw=vec(points[2],origin);
  const e2=unit(vec(raw,scale(e1,dot(raw,e1))));
  const e3=unit(cross(e1,e2));
  return {{o:origin,e1:e1,e2:e2,e3:e3}};
}}
function local(p,f) {{const v=vec(p,f.o); return {{x:dot(v,f.e1),y:dot(v,f.e2),z:dot(v,f.e3)}};}}
function world(v,f) {{return add(f.o,add(scale(f.e1,v.x),add(scale(f.e2,v.y),scale(f.e3,v.z))));}}
function overlayByThreeResidues() {{
  const qr=unique3(DATA.qres), hr=unique3(DATA.hres);
  if(!qr||!hr) throw new Error('Exactly three SPRITE residues are required for overlay');
  const qp=[],hp=[];
  for(let i=0;i<3;i++) {{
    const q=ca(qModel,DATA.qchain,qr[i]), h=ca(hModel,DATA.hchain,hr[i]);
    if(!q||!h) throw new Error('Could not locate a CA atom for all three SPRITE residues');
    qp.push({{x:q.x,y:q.y,z:q.z}}); hp.push({{x:h.x,y:h.y,z:h.z}});
  }}
  const qf=frame(qp), hf=frame(hp);
  for(const a of hModel.atoms) {{ const p=world(local({{x:a.x,y:a.y,z:a.z}},hf),qf); a.x=p.x;a.y=p.y;a.z=p.z; }}
  for(let i=0;i<3;i++) {{
    const h=ca(hModel,DATA.hchain,hr[i]), q=ca(qModel,DATA.qchain,qr[i]);
    if(!h||!q||Math.hypot(h.x-q.x,h.y-q.y,h.z-q.z)>0.15) throw new Error('SPRITE overlay verification failed');
  }}
}}
function styleFull() {{
  qModel.setStyle({{}},{{cartoon:{{color:'#b8bec7'}}}});
  hModel.setStyle({{}},{{cartoon:{{color:'#b8bec7'}}}});
}}
function styleOverlay() {{
  qModel.setStyle({{}},{{cartoon:{{color:BLUE}}}});
  hModel.setStyle({{}},{{cartoon:{{color:YELLOW}}}});
}}
function showSiteStyles() {{
  const qr=unique3(DATA.qres),hr=unique3(DATA.hres);
  if(!qr||!hr) throw new Error('Exactly three SPRITE residues are required');
  qModel.setStyle({{chain:DATA.qchain,resi:qr}},{{stick:{{radius:0.25,color:SITE_BLUE}},sphere:{{radius:0.42,color:SITE_BLUE}}}});
  hModel.setStyle({{chain:DATA.hchain,resi:hr}},{{stick:{{radius:0.25,color:SITE_YELLOW}},sphere:{{radius:0.42,color:SITE_YELLOW}}}});
}}
function showMode(mode) {{
  if(!viewer||!qModel||!hModel) return;
  try {{
    if(mode===3||mode===4) overlayByThreeResidues();
    if(mode===1||mode===2) styleFull(); else styleOverlay();
    if(mode===2||mode===4) {{ qModel.setStyle({{}},{{}}); hModel.setStyle({{}},{{}}); }}
    showSiteStyles();
    viewer.zoomTo(); viewer.render();
    const labels={{1:'Mode 1 — full proteins + 3 highlighted residues',2:'Mode 2 — 3 residues only',3:'Mode 3 — rigid overlay + highlighted residues',4:'Mode 4 — rigid overlay, 3 residues only'}};
    document.getElementById('status').textContent=labels[mode];
  }} catch(e) {{ document.getElementById('status').textContent='3D viewer error: '+e.message; }}
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
  window.applyMode=showMode;
  showMode(DATA.mode);
 }} catch(e) {{ document.getElementById('status').textContent='3D viewer error: '+e.message; }}
}}
init();
</script></body></html>'''
