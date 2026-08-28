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
        "query": query.upper(), "qchain": qchain, "qres": qres,
        "homolog": homolog.upper(), "hchain": hchain, "hres": hres,
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
let viewer=null,qModel=null,hModel=null;
const GREY='#a7adb5',BLUE='#2563eb',YELLOW='#facc15';
const SITE_BLUE='#2563eb',SITE_YELLOW='#f59e0b';

function unique3(arr) {{
 const out=[],seen=new Set();
 for(const x of (arr||[])) {{ const s=String(x); if(!seen.has(s)){{seen.add(s);out.push(s);}} }}
 return out.length===3?out:null;
}}
function ca(model,chain,res) {{
 const atoms=model.selectedAtoms({{chain:chain,resi:String(res)}});
 return atoms.find(a=>String(a.atom||'').toUpperCase()==='CA'||String(a.name||'').toUpperCase()==='CA')||null;
}}
function pt(a) {{ return {{x:Number(a.x),y:Number(a.y),z:Number(a.z)}}; }}
function sub(a,b) {{ return {{x:a.x-b.x,y:a.y-b.y,z:a.z-b.z}}; }}
function add(a,b) {{ return {{x:a.x+b.x,y:a.y+b.y,z:a.z+b.z}}; }}
function scale(a,s) {{ return {{x:a.x*s,y:a.y*s,z:a.z*s}}; }}
function dot(a,b) {{ return a.x*b.x+a.y*b.y+a.z*b.z; }}
function norm(a) {{ return Math.hypot(a.x,a.y,a.z); }}
function centroid(ps) {{ let c={{x:0,y:0,z:0}}; for(const p of ps)c=add(c,p); return scale(c,1/ps.length); }}

/* Horn's quaternion method: best-fit rigid rotation for the three matched CA atoms.
   Unlike the old three-point frame, this does not assume the two local triangles have
   identical side lengths. A rigid overlay can therefore be shown even when the
   experimental structures have locally different geometry. */
function bestFitRigid(src,dst) {{
 const cs=centroid(src),cd=centroid(dst);
 const P=src.map(p=>sub(p,cs)),Q=dst.map(p=>sub(p,cd));
 let Sxx=0,Sxy=0,Sxz=0,Syx=0,Syy=0,Syz=0,Szx=0,Szy=0,Szz=0;
 for(let i=0;i<3;i++) {{
   const p=P[i],q=Q[i];
   Sxx+=p.x*q.x; Sxy+=p.x*q.y; Sxz+=p.x*q.z;
   Syx+=p.y*q.x; Syy+=p.y*q.y; Syz+=p.y*q.z;
   Szx+=p.z*q.x; Szy+=p.z*q.y; Szz+=p.z*q.z;
 }}
 const N=[
  [Sxx+Syy+Szz,Syz-Szy,Szx-Sxz,Sxy-Syx],
  [Syz-Szy,Sxx-Syy-Szz,Sxy+Syx,Szx+Sxz],
  [Szx-Sxz,Sxy+Syx,-Sxx+Syy-Szz,Syz+Szy],
  [Sxy-Syx,Szx+Sxz,Syz+Szy,-Sxx-Syy+Szz]
 ];
 let q=[1,0,0,0];
 for(let k=0;k<80;k++) {{
   const n=[0,0,0,0];
   for(let r=0;r<4;r++)for(let c=0;c<4;c++)n[r]+=N[r][c]*q[c];
   const len=Math.hypot(n[0],n[1],n[2],n[3])||1;
   q=n.map(v=>v/len);
 }}
 const w=q[0],x=q[1],y=q[2],z=q[3];
 const R=[
  [1-2*(y*y+z*z),2*(x*y-w*z),2*(x*z+w*y)],
  [2*(x*y+w*z),1-2*(x*x+z*z),2*(y*z-w*x)],
  [2*(x*z-w*y),2*(y*z+w*x),1-2*(x*x+y*y)]
 ];
 function transform(p) {{
   const v=sub(p,cs);
   return add(cd,{{x:R[0][0]*v.x+R[0][1]*v.y+R[0][2]*v.z,y:R[1][0]*v.x+R[1][1]*v.y+R[1][2]*v.z,z:R[2][0]*v.x+R[2][1]*v.y+R[2][2]*v.z}});
 }}
 return {{transform:transform,src:src,dst:dst}};
}}

function getSpritePoints() {{
 const qr=unique3(DATA.qres),hr=unique3(DATA.hres);
 if(!qr||!hr)throw new Error('SPRITE requires exactly three residues on each structure');
 const qp=[],hp=[];
 for(let i=0;i<3;i++) {{
   const q=ca(qModel,DATA.qchain,qr[i]),h=ca(hModel,DATA.hchain,hr[i]);
   if(!q||!h)throw new Error('Could not locate CA atoms for all three SPRITE residues');
   qp.push(pt(q));hp.push(pt(h));
 }}
 return {{qr:qr,hr:hr,qp:qp,hp:hp}};
}}

function overlayHomolog() {{
 const {{qp,hp}}=getSpritePoints();
 const fit=bestFitRigid(hp,qp);
 const coords=hModel.atoms.map(a=>{{const p=fit.transform(pt(a));return [p.x,p.y,p.z];}});
 hModel.setCoordinates([coords],'array');
 const check=getSpritePoints();
 let sum=0,maxErr=0;
 for(let i=0;i<3;i++){{const d=Math.hypot(check.qp[i].x-check.hp[i].x,check.qp[i].y-check.hp[i].y,check.qp[i].z-check.hp[i].z);sum+=d*d;maxErr=Math.max(maxErr,d);}}
 const rmsd=Math.sqrt(sum/3);
 // Never hide a usable viewer merely because experimental local geometry differs.
 // The displayed value is diagnostic; the transformation is still the best-fit rigid fit.
 return {{rmsd:rmsd,maxErr:maxErr}};
}}

function baseGrey() {{
 qModel.setStyle({{}},{{cartoon:{{color:GREY}}}});
 hModel.setStyle({{}},{{cartoon:{{color:GREY}}}});
}}
function baseOverlay() {{
 qModel.setStyle({{}},{{cartoon:{{color:BLUE}}}});
 hModel.setStyle({{}},{{cartoon:{{color:YELLOW}}}});
}}
function highlightSites() {{
 const qr=unique3(DATA.qres),hr=unique3(DATA.hres);
 if(!qr||!hr)throw new Error('Exactly three SPRITE residues are required');
 qModel.addStyle({{chain:DATA.qchain,resi:qr}},{{stick:{{radius:0.30,color:SITE_BLUE}},sphere:{{radius:0.48,color:SITE_BLUE}}}});
 hModel.addStyle({{chain:DATA.hchain,resi:hr}},{{stick:{{radius:0.30,color:SITE_YELLOW}},sphere:{{radius:0.48,color:SITE_YELLOW}}}});
 // Explicit labels make it visually unambiguous that all three residues are present.
 qr.forEach((r,i)=>{{const a=ca(qModel,DATA.qchain,r);if(a)viewer.addLabel(DATA.qchain+':'+r,{{position:a,fontSize:11,fontColor:BLUE,backgroundColor:'white',backgroundOpacity:0.65}});}});
 hr.forEach((r,i)=>{{const a=ca(hModel,DATA.hchain,r);if(a)viewer.addLabel(DATA.hchain+':'+r,{{position:a,fontSize:11,fontColor:'#9a6700',backgroundColor:'white',backgroundOpacity:0.65}});}});
}}
function siteOnly() {{
 const qr=unique3(DATA.qres),hr=unique3(DATA.hres);
 qModel.setStyle({{}},{{}});hModel.setStyle({{}},{{}});
 qModel.setStyle({{chain:DATA.qchain,resi:qr}},{{stick:{{radius:0.32,color:SITE_BLUE}},sphere:{{radius:0.50,color:SITE_BLUE}}}});
 hModel.setStyle({{chain:DATA.hchain,resi:hr}},{{stick:{{radius:0.32,color:SITE_YELLOW}},sphere:{{radius:0.50,color:SITE_YELLOW}}}});
}}
function clearLabels() {{ viewer.removeAllLabels(); }}
function showMode(mode) {{
 if(!viewer||!qModel||!hModel)return;
 try {{
   clearLabels();
   let diag=null;
   if(mode===3||mode===4)diag=overlayHomolog();
   if(mode===1){{baseGrey();highlightSites();}}
   else if(mode===2){{siteOnly();}}
   else if(mode===3){{baseOverlay();highlightSites();}}
   else {{siteOnly();}}
   viewer.zoomTo();viewer.render();
   const labels={{1:'Mode 1 — full proteins grey + three sites highlighted',2:'Mode 2 — three SPRITE residues only',3:'Mode 3 — best-fit rigid overlay + three sites highlighted',4:'Mode 4 — best-fit rigid overlay, three sites only'}};
   let text=labels[mode];
   if(diag)text+=' · 3-site CA RMSD '+diag.rmsd.toFixed(2)+' Å';
   document.getElementById('status').textContent=text;
 }}catch(e){{document.getElementById('status').textContent='3D viewer error: '+e.message;}}
}}
async function loadPdb(id) {{
 const r=await fetch('https://files.rcsb.org/download/'+encodeURIComponent(id)+'.pdb');
 if(!r.ok)throw new Error('Could not download '+id+' from RCSB');
 return await r.text();
}}
async function init() {{
 try {{
   viewer=$3Dmol.createViewer(document.getElementById('viewer'),{{backgroundColor:'white'}});
   const qp=await loadPdb(DATA.query);qModel=viewer.addModel(qp,'pdb');
   if(!DATA.homolog){{qModel.setStyle({{}},{{cartoon:{{color:GREY}}}});viewer.zoomTo();viewer.render();document.getElementById('status').textContent='Select a SPRITE homolog to compare structures.';return;}}
   const hp=await loadPdb(DATA.homolog);hModel=viewer.addModel(hp,'pdb');showMode(DATA.mode);
 }}catch(e){{document.getElementById('status').textContent='3D viewer error: '+e.message;}}
}}
init();
</script></body></html>'''
