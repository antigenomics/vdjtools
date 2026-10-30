"""Build the README's schematic analysis map using deterministic SVG primitives.

Palette follows isalgo.dev. Diagrams illustrate analysis types; motif inputs are a recorded public VDJdb example.
Render PNG with: rsvg-convert assets/analysis-atlas-light.svg -o assets/analysis-atlas.png
"""
from pathlib import Path
from html import escape
import csv
import io
import math
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
for dark in (False, True):
    theme = 'dark' if dark else 'light'
    bg, ink, soft, line = ('#1a2332','#eef2f7','#9cadbf','#354254') if dark else ('#ffffff','#1a2332','#54657a','#e6e9ee')
    teal = '#25c39d' if dark else '#16a085'
    amber = '#ffbd59' if dark else '#f0a93f'
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="1430" viewBox="0 0 1200 1430" role="img" aria-labelledby="title desc"><title id="title">VDJtools: immune-repertoire analysis</title><desc id="desc">A schematic workflow from bulk and paired receptor tables to nine analysis types. Charts are schematic except the recorded VDJdb motif example; the logo shows observed-residue frequency differences against V/J/length-matched background.</desc><rect width="1200" height="1430" fill="{bg}"/><g font-family="Arial, Helvetica, sans-serif">']
    def text(x,y,s,size=20,color=None,anchor='start',weight='400'):
        out.append(f'<text x="{x}" y="{y}" font-size="{size}" fill="{color or ink}" text-anchor="{anchor}" font-weight="{weight}">{escape(s)}</text>')
    def path(d,c=None,w=2,opacity=1):
        out.append(f'<path d="{d}" fill="none" stroke="{c or soft}" stroke-width="{w}" opacity="{opacity}" stroke-linecap="round" stroke-linejoin="round"/>')
    def rect(x,y,w,h,c,opacity=1):
        out.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="{c}" opacity="{opacity}"/>')
    def circle(x,y,r,c,fill=None):
        out.append(f'<circle cx="{x}" cy="{y}" r="{r}" fill="{fill or bg}" stroke="{c}" stroke-width="2"/>')
    text(65,55,'ISALGO  /  VDJTOOLS',16,teal,weight='700')
    text(65,108,'Immune-repertoire analysis',30,weight='700')
    text(65,145,'From receptor tables to sequence, sample and cohort measurements',18,soft)
    # Original vector logo, embedded without rasterization or external references.
    logo=ET.parse(ROOT/f'assets/vdjtools_{theme}.svg').getroot()
    logo.set('x','925');logo.set('y','18');logo.set('width','230');logo.set('height','145')
    out.append(ET.tostring(logo,encoding='unicode'))
    path('M65 176 H1135',line)
    # Inputs: sequence rows and paired-chain barcodes, converging on a shared table.
    for i in range(4):
        y=219+i*16
        for x,w,c in [(75,65,soft),(145,15,amber),(165,50,teal)]:rect(x,y,w,5,c)
    text(250,233,'Bulk repertoires',23,weight='700');text(250,266,'AIRR · MIGEC / ARDA · immunoSEQ',15,soft)
    for i in range(3):
        x=695+i*36;circle(x,244,17,soft)
        path(f'M{x-8} 237 l8 9 8 -9 M{x} 246 v9',teal,2)
    text(810,233,'Paired receptors',23,weight='700');text(810,266,'10x / CellRanger · AIRR cells',15,soft)
    path('M400 290 V311 H600 M935 290 V311 H600 V338',soft)
    text(600,370,'Read  →  validate  →  filter  →  normalize',23,anchor='middle',weight='700')
    text(600,401,'Clonotype tables or cell-aware chain tables',18,soft,anchor='middle')
    path('M600 419 V438 M230 438 H990 M230 438 V457 M610 438 V457 M990 438 V457',line,2)
    titles=['Usage & spectratypes','Diversity & sampling','Sequence motifs','Recombination models','Overlap & enrichment','Normalization & bias','Paired-cell analysis','Cohort representations','Downstream tools']
    subs=[('V / D / J / C usage · VJ pairs','Junction lengths · isotype summaries'),('Clone sizes · richness · diversity','Rarefaction · coverage'),('k-mer cloud · corrected PWM','Hydropathy · physicochemistry · SHM'),('Pgen · sequence generation','Latent-history model fitting'),('Repertoire overlap · distances','TCRnet · biomarkers · dynamics'),('Depth normalization · filtering','Cross-batch V/J-usage correction'),('Chain QC · pairing · paired Pgen','Single-cell interoperability'),('Statistical features · sample PCA','Corpus signatures · vsig / rsig'),('VDJmatch · VDJdb queries','mirpy · receptor & repertoire embeddings')]
    for k in range(9):
        x=45+(k%3)*380;y=480+(k//3)*300
        out.append(f'<g transform="translate({x},{y})">')
        text(0,0,f'{k+1:02}',15,teal,weight='700')
        if k==0:
            path('M25 135 H310 M25 135 V30',line)
            for j,h in enumerate([22,40,75,95,71,45,27,15]):rect(40+j*23,135-h,14,h,teal,0.35+j*.07)
            path('M38 126 C70 125 84 90 115 63 S156 56 181 89 S224 124 260 126',soft,2.5)
        elif k==1:
            path('M25 30 V135 H310',soft)
            for d,c in [('M25 135 C70 65 170 36 302 33',teal),('M25 135 C65 83 175 62 302 58',soft),('M25 135 C70 100 178 92 302 89',amber)]:path(d,c,3)
        elif k==2:
            # Recorded public motif frequencies, not decorative random cells.
            with (ROOT/'assets/source/kmer-example.tsv').open() as f:
                words = list(csv.DictReader(f, delimiter='\t'))
            for row, (cx, cy) in zip(words, [(65,40),(55,67),(72,94),(21,17),(107,17),(18,112),(113,112)]):
                size = 12 + 14 * int(row['count']) / max(int(r['count']) for r in words)
                text(cx, cy, row['kmer'], size, teal if int(row['count']) > 1 else soft,
                     anchor='middle', weight='700')
            # Signed frequency-difference logo: p_motif(aa, position) - p_background.
            # The source table's background is matched on V, J and junction length.
            with (ROOT/'assets/source/motif-example.tsv').open() as f:
                pwm = list(csv.DictReader(f, delimiter='\t'))
            # Residue classes share a color on both sides of the zero baseline.
            residue_colors = {
                aa: (night if dark else day)
                for residues, day, night in [
                    ('G', '#b17b20', '#e5b65b'),
                    ('P', '#946844', '#cba784'),
                    ('C', '#a48a21', '#d9c45f'),
                    ('ST', '#168b78', '#43c4ad'),
                    ('NQ', '#8063a6', '#b99adf'),
                    ('DE', '#c05c4b', '#ef9586'),
                    ('KRH', '#3f75ac', '#83b5e7'),
                    ('FWY', '#98658c', '#d3a0c7'),
                    ('AILMV', '#54657a', '#b5c3d3'),
                ] for aa in residues
            }
            path('M153 75 H335', line)
            for pos in range(13):
                column = [(r['aa'], float(r['freq']) - float(r['freq.bg']))
                          for r in pwm if int(r['pos']) == pos]
                for sign in (1, -1):
                    offset = 0.0
                    for aa, delta in sorted(column, key=lambda v: abs(v[1])):
                        if delta * sign <= 0:
                            continue
                        height = abs(delta) * 50
                        baseline = 75 - offset if sign > 0 else 75 + offset + height
                        # Arial capital height is 73 units at font-size 100.
                        out.append(f'<text transform="translate({155+pos*14},{baseline}) scale(0.17,{height/73})" '
                                   f'font-size="100" '
                                   f'fill="{residue_colors[aa]}" font-weight="700">{aa}</text>')
                        offset += height
            text(244,108,'Δ frequency vs V/J background',12,soft,anchor='middle')
            # Use the shipped hydropathy scale, exactly as physchem_profile does.
            raw = (ROOT/'python/vdjtools/resources/aa_property_table.txt').read_text()
            table = '\n'.join(ln for ln in raw.splitlines() if ln and not ln.startswith('##'))
            hydropathy = {r['amino_acid']: float(r['hydropathy'])
                          for r in csv.DictReader(io.StringIO(table), delimiter='\t')}
            consensus = [max((r for r in pwm if int(r['pos']) == pos),
                             key=lambda r: float(r['freq']))['aa'] for pos in range(13)]
            for pos, aa in enumerate(consensus):
                fraction = (hydropathy[aa] + 4.5) / 9
                rgb = tuple(round(a + fraction * (b-a))
                            for a, b in zip((100,170,165), (236,180,92)))
                rect(22+pos*24,120,22,20,f'rgb{rgb}')
                text(33+pos*24,135,aa,17,'#1a2332',anchor='middle',weight='700')
            text(178,155,'Hydropathy: hydrophilic → hydrophobic',11,soft,anchor='middle')
        elif k==3:
            for j,(s,c) in enumerate([('V',soft),('D',amber),('J',teal)]):
                rect(27+j*95,30,65,16,c);text(60+j*95,23,s,16,c,anchor='middle')
                path(f'M{60+j*95} 52 V68 L155 91',c)
            for j,(w,c) in enumerate([(85,soft),(22,amber),(65,teal)]):
                start=[60,149,175][j];rect(start,96,w,12,c)
            text(155,140,'V(D)J recombination',17,soft,anchor='middle')
        elif k==4:
            pts=[(55,65),(88,35),(115,83),(67,117),(146,49),(180,94),(220,37),(254,84),(225,126),(288,49)]
            for a,b in [(0,1),(0,2),(0,3),(1,2),(1,4),(2,3),(2,4),(2,5),(4,5),(5,6),(5,7),(6,7),(6,9),(7,8),(7,9)]:
                path(f'M{pts[a][0]} {pts[a][1]} L{pts[b][0]} {pts[b][1]}',line,2)
            for j,(a,b) in enumerate(pts):circle(a,b,6,teal if j<6 else soft,teal if j<6 else bg)
        elif k==5:
            for start,hs in [(28,[25,80,45,100]),(198,[55,65,58,70])]:
                path(f'M{start-5} 135 H{start+105}',soft)
                for j,h in enumerate(hs):rect(start+j*25,135-h,16,h,teal if j%2 else soft,.7)
            path('M137 78 H179 M172 72 L179 78 172 84',soft)
        elif k==6:
            for cx in [82,240]:
                circle(cx,83,48,line);circle(cx,83,22,line)
                path(f'M{cx-24} 54 v-20 m0 12 l-9 -9 m9 9 l9 -9',teal,3)
                path(f'M{cx+24} 54 v-20 m0 12 l-9 -9 m9 9 l9 -9',soft,3)
            path('M113 35 C150 8 173 8 210 35',teal)
            for j in range(17):rect(126+j*4,100,1+(j%2),24,soft)
        elif k==7:
            for i in range(5):
                for j in range(5):rect(24+j*17,40+i*17,14,14,teal,((i*3+j*7)%9+1)/10)
            path('M126 83 H154 M147 77 L154 83 147 89',soft)
            path('M175 30 V135 H306',line)
            for j in range(24):
                cx=202+(j%2)*65+math.sin(j*4)*17;cy=97-(j%2)*45+math.cos(j*3)*17
                circle(round(cx,2),round(cy,2),3,teal if j%2 else soft,teal if j%2 else soft)
        else:
            for j in range(4):
                rect(15,60+j*12,48,4,soft)
            path('M77 80 H100 V42 H124 M100 80 V119 H124',soft)
            text(138,48,'VDJmatch',21,teal,weight='700')
            text(138,125,'mirpy',21,teal,weight='700')
            path('M253 42 H274 M268 37 L274 42 268 47',soft)
            out.append(f'<ellipse cx="302" cy="30" rx="20" ry="7" fill="{bg}" stroke="{soft}" stroke-width="2"/>')
            path('M282 30 V58 C282 68 322 68 322 58 V30 M282 45 C282 55 322 55 322 45',soft)
            path('M219 119 H250 M244 114 L250 119 244 124',soft)
            for j, h in enumerate([18,30,12,26,21]):
                rect(271+j*11,136-h,6,h,teal if j%2 else soft)
        text(0,180,titles[k],23,weight='700')
        text(0,211,subs[k][0],16,soft);text(0,237,subs[k][1],16,soft)
        out.append('</g>')
    path('M45 1360 H1155',line)
    text(600,1400,'Python library  ·  C++ acceleration  ·  CLI  ·  multiprocessing  ·  Polars  ·  Parquet',18,soft,anchor='middle')
    out.append('</g></svg>')
    (ROOT/f'assets/analysis-atlas-{theme}.svg').write_text('\n'.join(out)+'\n')
