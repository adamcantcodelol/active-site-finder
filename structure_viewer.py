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
    q_aln: str | None = None,
    t_aln: str | None = None,
    q_start: int | None = None,
    t_start: int | None = None,
) -> str:
    """Create a self-contained 3Dmol viewer with four SPRITE comparison modes.

    Modes 3/4 use the complete residue correspondence supplied by Foldseek's
    structural alignment, rather than trying to orient two proteins from only
    the three SPRITE residues. The three SPRITE residues are then highlighted
    on top of that scientifically meaningful structural overlay.
    """
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
    if homolog and (not q_aln or not t_aln):
        raise ValueError("The structural homolog is missing Foldseek's full alignment")
    if mode not in (1, 2, 3, 4):
        mode = 1

    safe_height = max(420, min(int(height), 900))
    payload = {
        "query": query.upper(), "qchain": qchain, "qres": qres,
        "homolog": homolog.upper(), "hchain": hchain, "hres": hres,
        "mode": mode, "q_aln": q_aln or "", "t_aln": t_aln or "",
        "q_start": q_start, "t_start": t_start,
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
function centroid(ps) {{ let c={{x:0,y:0,z:0}}; for(const p of ps)c=add(c,p); return scale(c,1/ps.length); }}

function caResidues(model,chain) {{
 const atoms=model.selectedAtoms({{chain:chain}});
 const out=[],seen=new Set();
 for(const a of atoms) {{
   const name=String(a.atom||a.name||'').toUpperCase();
   if(name!=='CA')continue;
   const key=String(a.resi)+(a.resn?':'+a.resn:'');
   if(seen.has(key))continue;
   seen.add(key);out.push({{resnum:String(a.resi),resn:String(a.resn||''),atom:a}});
 }}
 return out;
}}

function alignmentStart(aln,residues,startHint) {{
 const probe=String(aln||'').split('').filter(c=>c!=='-').slice(0,80).map(c=>c.toUpperCase());
 if(!probe.length||!residues.length)return -1;
 const hint=Math.max(0,(Number(startHint)||1)-1);
 let best=-1,bestScore=-1e99;
 for(let i=0;i<residues.length;i++) {{
   let matches=0;
   for(let j=0;j<probe.length && i+j<residues.length;j++) {{
     const aa=oneLetter(residues[i+j].resn);
     if(aa===probe[j]||probe[j]==='X')matches++;
   }}
   const score=matches-Math.min(Math.abs(i-hint),200)*0.002;
   if(score>bestScore){{bestScore=score;best=i;}}
 }}
 const identity=best>=0?bestScore/Math.max(1,probe.length):0;
 return identity>=0.60?best:-1;
}}
function oneLetter(resn) {{
 const m={{ALA:'A',ARG:'R',ASN:'N',ASP:'D',CYS:'C',GLN:'Q',GLU:'E',GLY:'G',HIS:'H',ILE:'I',LEU:'L',LYS:'K',MET:'M',PHE:'F',PRO:'P',SER:'S',THR:'T',TRP:'W',TYR:'Y',VAL:'V',MSE:'M',SEP:'S',TPO:'T',PTR:'Y'}};
 return m[String(resn||'').toUpperCase()]||'X';
}}

/* Build residue-to-residue correspondence from Foldseek's complete alignment.
   This is the same alignment Foldseek used to calculate structural similarity;
   we are not inventing a new sequence alignment in the browser. */
function getFullAlignmentPairs() {{
 const qr=caResidues(qModel,DATA.qchain),hr=caResidues(hModel,DATA.hchain);
 let qi=alignmentStart(DATA.q_aln,qr,DATA.q_start),hi=alignmentStart(DATA.t_aln,hr,DATA.t_start);
 if(qi<0||hi<0)throw new Error('Could not anchor Foldseek alignment to deposited C-alpha coordinates');
 const pairs=[];
 for(let k=0;k<Math.min(DATA.q_aln.length,DATA.t_aln.length);k++) {{
   const qc=DATA.q_aln[k],tc=DATA.t_aln[k];
   if(qc!=='-'&&tc!=='-'&&qi<qr.length&&hi<hr.length) {{
     const q=qr[qi],h=hr[hi];
     if((qc.toUpperCase()==='X'||oneLetter(q.resn)===qc.toUpperCase())&&(tc.toUpperCase()==='X'||oneLetter(h.resn)===tc.toUpperCase())) {{
       pairs.push({{q:q.atom,h:h.atom,qres:q.resnum,hres:h.resnum}});
     }}
   }}
   if(qc!=='-')qi++;
   if(tc!=='-')hi++;
 }}
 if(pairs.length<3)throw new Error('Foldseek alignment contains fewer than three usable C-alpha pairs');
 return pairs;
}}

/* Horn quaternion best-fit rigid transform over the FULL Foldseek alignment.
   This performs only rotation + translation: no scaling and no deformation. */
function bestFitRigid(src,dst) {{
 const cs=centroid(src),cd=centroid(dst);
 const P=src.map(p=>sub(p,cs)),Q=dst.map(p=>sub(p,cd));
 let Sxx=0,Sxy=0,Sxz=0,Syx=0,Syy=0,Syz=0,Szx=0,Szy=0,Szz=0;
 for(let i=0;i<P.length;i++) {{ const p=P[i],q=Q[i];
   Sxx+=p.x*q.x;Sxy+=p.x*q.y;Sxz+=p.x*q.z;
   Syx+=p.y*q.x;Syy+=p.y*q.y;Syz+=p.y*q.z;
   Szx+=p.z*q.x;Szy+=p.z*q.y;Szz+=p.z*q.z; }}
 const N=[[Sxx+Syy+Szz,Syz-Szy,Szx-Sxz,Sxy-Syx],[Syz-Szy,Sxx-Syy-Szz,Sxy+Syx,Szx+Sxz],[Szx-Sxz,Sxy+Syx,-Sxx+Syy-Szz,Syz+Szy],[Sxy-Syx,Szx+Sxz,Syz+Szy,-Sxx-Syy+Szz]];
 let v=[1,0,0,0];
 for(let k=0;k<100;k++){{let n=[0,0,0,0];for(let r=0;r<4;r++)for(let c=0;c<4;c++)n[r]+=N[r][c]*v[c];const l=Math.hypot(...n)||1;v=n.map(x=>x/l);}}
 const w=v[0],x=v[1],y=v[2],z=v[3];
 const R=[[1-2*(y*y+z*z),2*(x*y-w*z),2*(x*z+w*y)],[2*(x*y+w*z),1-2*(x*x+z*z),2*(y*z-w*x)],[2*(x*z-w*y),2*(y*z+w*x),1-2*(x*x+y*y)]];
 return function(p) {{const q=sub(p,cs);return add(cd,{{x:R[0][0]*q.x+R[0][1]*q.y+R[0][2]*q.z,y:R[1][0]*q.x+R[1][1]*q.y+R[1][2]*q.z,z:R[2][0]*q.x+R[2][1]*q.y+R[2][2]*q.z}});}};
}}

function overlayByFoldseek() {{
 const pairs=getFullAlignmentPairs();
 const src=pairs.map(p=>pt(p.h)),dst=pairs.map(p=>pt(p.q));
 const transform=bestFitRigid(src,dst);
 const coords=hModel.atoms.map(a=>{{const p=transform(pt(a));return [p.x,p.y,p.z];}});
 hModel.setCoordinates([coords],'array');
 let sum=0;
 for(let i=0;i<pairs.length;i++){{const p=transform(src[i]),q=dst[i];sum+=(p.x-q.x)**2+(p.y-q.y)**2+(p.z-q.z)**2;}}
 return {{rmsd:Math.sqrt(sum/pairs.length),pairs:pairs}};
}}

function baseGrey() {{ qModel.setStyle({{}},{{cartoon:{{color:GREY}}}});hModel.setStyle({{}},{{cartoon:{{color:GREY}}}}); }}
function baseOverlay() {{ qModel.setStyle({{}},{{cartoon:{{color:BLUE}}}});hModel.setStyle({{}},{{cartoon:{{color:YELLOW}}}}); }}
function highlightSites() {{
 const qr=unique3(DATA.qres),hr=unique3(DATA.hres);if(!qr||!hr)throw new Error('Exactly three SPRITE residues are required');
 qModel.addStyle({{chain:DATA.qchain,resi:qr}},{{stick:{{radius:0.30,color:SITE_BLUE}},sphere:{{radius:0.48,color:SITE_BLUE}}}});
 hModel.addStyle({{chain:DATA.hchain,resi:hr}},{{stick:{{radius:0.30,color:SITE_YELLOW}},sphere:{{radius:0.48,color:SITE_YELLOW}}}});
}}
function siteOnly() {{
 const qr=unique3(DATA.qres),hr=unique3(DATA.hres);if(!qr||!hr)throw new Error('Exactly three SPRITE residues are required');
 qModel.setStyle({{}},{{}});hModel.setStyle({{}},{{}});
 qModel.setStyle({{chain:DATA.qchain,resi:qr}},{{stick:{{radius:0.34,color:SITE_BLUE}},sphere:{{radius:0.52,color:SITE_BLUE}}}});
 hModel.setStyle({{chain:DATA.hchain,resi:hr}},{{stick:{{radius:0.34,color:SITE_YELLOW}},sphere:{{radius:0.52,color:SITE_YELLOW}}}});
}}
function clearLabels() {{ viewer.removeAllLabels(); }}
function addSiteLabels() {{
 const qr=unique3(DATA.qres),hr=unique3(DATA.hres);
 qr.forEach(r=>{{const a=ca(qModel,DATA.qchain,r);if(a)viewer.addLabel(DATA.qchain+':'+r,{{position:a,fontSize:11,fontColor:BLUE,backgroundColor:'white',backgroundOpacity:0.65}});}});
 hr.forEach(r=>{{const a=ca(hModel,DATA.hchain,r);if(a)viewer.addLabel(DATA.hchain+':'+r,{{position:a,fontSize:11,fontColor:'#9a6700',backgroundColor:'white',backgroundOpacity:0.65}});}});
}}
function showMode(mode) {{
 if(!viewer||!qModel||!hModel)return;
 try {{
   clearLabels();let diag=null;
   if(mode===3||mode===4)diag=overlayByFoldseek();
   if(mode===1){{baseGrey();highlightSites();}}
   else if(mode===2){{siteOnly();}}
   else if(mode===3){{baseOverlay();highlightSites();addSiteLabels();}}
   else {{siteOnly();addSiteLabels();}}
   viewer.zoomTo();viewer.render();
   const labels={{1:'Mode 1 — full proteins grey + three sites highlighted',2:'Mode 2 — three SPRITE residues only',3:'Mode 3 — Foldseek-aligned full overlay + three sites highlighted',4:'Mode 4 — Foldseek-aligned overlay, three sites only'}};
   let text=labels[mode];if(diag)text+=' · full-alignment CA RMSD '+diag.rmsd.toFixed(2)+' Å · '+diag.pairs.length+' aligned residues';
   document.getElementById('status').textContent=text;
 }}catch(e){{document.getElementById('status').textContent='3D viewer error: '+e.message;}}
}}
async function loadPdb(id) {{const r=await fetch('https://files.rcsb.org/download/'+encodeURIComponent(id)+'.pdb');if(!r.ok)throw new Error('Could not download '+id+' from RCSB');return await r.text();}}
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
