"""Interactive 3D structure comparison helpers for MBRC Active Site Finder."""
from __future__ import annotations

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
    payload = {
        "query": query.upper(),
        "qchain": qchain,
        "qres": qres,
        "homolog": homolog.upper(),
        "hchain": hchain,
        "hres": hres,
        "mode": mode,
    }
    data = json.dumps(payload).replace("</", "<\\/")

    return f'''<!doctype html>
<html><head><meta charset="utf-8">
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>html,body{{margin:0;padding:0;background:#fff;font-family:Arial,Helvetica,sans-serif}}#viewer{{width:100%;height:{safe_height-52}px;position:relative}}#status{{height:52px;display:flex;align-items:center;padding:0 12px;box-sizing:border-box;color:#64748b;font-size:13px;border-top:1px solid #e2e8f0}}</style>
</head><body><div id="viewer"></div><div id="status">Loading structures…</div>
<script>
const DATA={data};
let viewer=null, qModel=null, hModel=null;
const GREY='#a7adb5', BLUE='#2563eb', YELLOW='#facc15';
const SITE_BLUE='#2563eb', SITE_YELLOW='#f59e0b';

function unique3(arr) {{
  const out=[]; const seen=new Set();
  for(const x of (arr||[])) {{
    const s=String(x);
    if(!seen.has(s)) {{ seen.add(s); out.push(s); }}
  }}
  return out.length===3 ? out : null;
}}

function ca(model, chain, res) {{
  const atoms=model.selectedAtoms({{chain:chain,resi:String(res)}});
  return atoms.find(a=>String(a.atom||'').toUpperCase()==='CA' || String(a.name||'').toUpperCase()==='CA') || null;
}}

function vsub(a,b) {{ return {{x:a.x-b.x,y:a.y-b.y,z:a.z-b.z}}; }}
function vadd(a,b) {{ return {{x:a.x+b.x,y:a.y+b.y,z:a.z+b.z}}; }}
function vscale(a,s) {{ return {{x:a.x*s,y:a.y*s,z:a.z*s}}; }}
function dot(a,b) {{ return a.x*b.x+a.y*b.y+a.z*b.z; }}
function cross(a,b) {{ return {{x:a.y*b.z-a.z*b.y,y:a.z*b.x-a.x*b.z,z:a.x*b.y-a.y*b.x}}; }}
function norm(a) {{ return Math.hypot(a.x,a.y,a.z); }}
function unit(a) {{ const n=norm(a); if(n<1e-8) throw new Error('The three SPRITE residues are geometrically degenerate'); return vscale(a,1/n); }}

/*
 * Build a right-handed orthonormal frame from the three CA atoms.
 * Mapping homolog frame -> query frame gives a rigid-body transformation:
 * translation + rotation only. No scaling, stretching, or camera tricks.
 */
function makeFrame(points) {{
  const o=points[0];
  const e1=unit(vsub(points[1],o));
  const raw=vsub(points[2],o);
  const e2raw=vsub(raw,vscale(e1,dot(raw,e1)));
  const e2=unit(e2raw);
  const e3=unit(cross(e1,e2));
  return {{o:o,e1:e1,e2:e2,e3:e3}};
}}

function toLocal(p,f) {{
  const v=vsub(p,f.o);
  return {{x:dot(v,f.e1),y:dot(v,f.e2),z:dot(v,f.e3)}};
}}

function toWorld(p,f) {{
  return vadd(f.o,vadd(vscale(f.e1,p.x),vadd(vscale(f.e2,p.y),vscale(f.e3,p.z))));
}}

function getSpritePoints() {{
  const qr=unique3(DATA.qres), hr=unique3(DATA.hres);
  if(!qr||!hr) throw new Error('SPRITE requires exactly three residues on each structure');
  const qp=[],hp=[];
  for(let i=0;i<3;i++) {{
    const q=ca(qModel,DATA.qchain,qr[i]);
    const h=ca(hModel,DATA.hchain,hr[i]);
    if(!q||!h) throw new Error('Could not locate CA atoms for all three SPRITE residues');
    qp.push({{x:Number(q.x),y:Number(q.y),z:Number(q.z)}});
    hp.push({{x:Number(h.x),y:Number(h.y),z:Number(h.z)}});
  }}
  return {{qr,hr,qp,hp}};
}}

function rigidOverlayThreeResidues() {{
  const {{qp,hp}}=getSpritePoints();
  const qf=makeFrame(qp), hf=makeFrame(hp);

  // Transform every homolog atom using the same rigid transform used for the site.
  const coords=hModel.atoms.map(a=>{{
    const local=toLocal({{x:Number(a.x),y:Number(a.y),z:Number(a.z)}},hf);
    const world=toWorld(local,qf);
    return [world.x,world.y,world.z];
  }});
  hModel.setCoordinates([coords],'array');
  hModel.molObj=null;

  // Verify the actual model coordinates after setCoordinates(), not a stale atom copy.
  const check=getSpritePoints();
  let maxErr=0;
  for(let i=0;i<3;i++) {{
    const d=Math.hypot(check.qp[i].x-check.hp[i].x,check.qp[i].y-check.hp[i].y,check.qp[i].z-check.hp[i].z);
    maxErr=Math.max(maxErr,d);
  }}
  if(maxErr>0.35) throw new Error('SPRITE overlay could not be aligned from all three matched residues (max CA error '+maxErr.toFixed(2)+' Å)');
  return maxErr;
}}

function styleFullGrey() {{
  qModel.setStyle({{}},{{cartoon:{{color:GREY}}}});
  hModel.setStyle({{}},{{cartoon:{{color:GREY}}}});
}}

function styleFullOverlay() {{
  qModel.setStyle({{}},{{cartoon:{{color:BLUE}}}});
  hModel.setStyle({{}},{{cartoon:{{color:YELLOW}}}});
}}

function highlightSites() {{
  const qr=unique3(DATA.qres), hr=unique3(DATA.hres);
  if(!qr||!hr) throw new Error('Exactly three SPRITE residues are required');

  // Add site styles instead of replacing the base cartoon. This makes all three
  // residues remain visible as part of the protein while being unmistakable.
  qModel.addStyle({{chain:DATA.qchain,resi:qr}},{{stick:{{radius:0.25,color:SITE_BLUE}},sphere:{{radius:0.38,color:SITE_BLUE}}}});
  hModel.addStyle({{chain:DATA.hchain,resi:hr}},{{stick:{{radius:0.25,color:SITE_YELLOW}},sphere:{{radius:0.38,color:SITE_YELLOW}}}});
}}

function showSiteOnly() {{
  const qr=unique3(DATA.qres), hr=unique3(DATA.hres);
  qModel.setStyle({{}},{{}});
  hModel.setStyle({{}},{{}});
  qModel.setStyle({{chain:DATA.qchain,resi:qr}},{{stick:{{radius:0.28,color:SITE_BLUE}},sphere:{{radius:0.46,color:SITE_BLUE}}}});
  hModel.setStyle({{chain:DATA.hchain,resi:hr}},{{stick:{{radius:0.28,color:SITE_YELLOW}},sphere:{{radius:0.46,color:SITE_YELLOW}}}});
}}

function showMode(mode) {{
  if(!viewer||!qModel||!hModel) return;
  try {{
    // Every iframe starts from the original structures, so there is no stale
    // transformed homolog when the Streamlit widget is rerun.
    if(mode===3||mode===4) rigidOverlayThreeResidues();

    if(mode===1) {{
      styleFullGrey();
      highlightSites();
    }} else if(mode===2) {{
      showSiteOnly();
    }} else if(mode===3) {{
      styleFullOverlay();
      highlightSites();
    }} else {{
      showSiteOnly();
    }}

    viewer.zoomTo();
    viewer.render();
    const labels={{
      1:'Mode 1 — full proteins grey + three sites highlighted',
      2:'Mode 2 — three SPRITE residues only',
      3:'Mode 3 — rigid three-residue overlay + sites highlighted',
      4:'Mode 4 — rigid three-residue overlay, sites only'
    }};
    document.getElementById('status').textContent=labels[mode];
  }} catch(e) {{
    document.getElementById('status').textContent='3D viewer error: '+e.message;
  }}
}}

async function loadPdb(id) {{
  const r=await fetch('https://files.rcsb.org/download/'+encodeURIComponent(id)+'.pdb');
  if(!r.ok) throw new Error('Could not download '+id+' from RCSB');
  return await r.text();
}}

async function init() {{
  try {{
    viewer=$3Dmol.createViewer(document.getElementById('viewer'),{{backgroundColor:'white'}});
    const qp=await loadPdb(DATA.query);
    qModel=viewer.addModel(qp,'pdb');
    if(!DATA.homolog) {{
      qModel.setStyle({{}},{{cartoon:{{color:GREY}}}});
      viewer.zoomTo(); viewer.render();
      document.getElementById('status').textContent='Select a SPRITE homolog to compare structures.';
      return;
    }}
    const hp=await loadPdb(DATA.homolog);
    hModel=viewer.addModel(hp,'pdb');
    showMode(DATA.mode);
  }} catch(e) {{
    document.getElementById('status').textContent='3D viewer error: '+e.message;
  }}
}}
init();
</script></body></html>'''
