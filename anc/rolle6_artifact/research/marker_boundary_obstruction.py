"""Further exact marker certificates, kept separate from the earlier H6 stage."""
import json
from pathlib import Path
from inflection_obstruction import patterns
from marker_hermite6_obstruction import check_one
from hermite_fiber_obstruction import check_full_identities
ROOT=Path(__file__).resolve().parents[1]

def check_certificates():
 check_full_identities();data=json.loads((ROOT/'data/marker_boundary_certificates.json').read_text());seen=set();leaves=0
 for c in data['certificates']:
  sig=tuple(c['signature']);assert sig not in seen;seen.add(sig);leaves+=check_one(c)['leaves']
 print('EXACTLY VERIFIED',len(seen),'further marker signatures and',leaves,'rational Bernstein rectangles.',flush=True)
 return seen

_masks=None
def reasons(word):
 global _masks
 if _masks is None:_masks=[c['mask']for c in json.loads((ROOT/'data/marker_boundary_certificates.json').read_text())['certificates']]
 sig=patterns(word)[0][0]
 return ['marker_boundary_elimination']if any(all(sig[i]==v for i,v in mask)for mask in _masks)else[]

if __name__=='__main__':check_certificates()
