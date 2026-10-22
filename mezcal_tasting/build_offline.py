"""Builds standalone/mezcal_tasting_offline.html with all data embedded. Needs server.js running on :3055."""
import json,urllib.request,urllib.parse,concurrent.futures as cf,os
P=json.load(open('data/palenques.json'))
def f(p):
    for _ in range(2):
        try:
            r=json.load(urllib.request.urlopen("http://localhost:3055/api/mezcals?name="+urllib.parse.quote(p['name'])+"&brand="+urllib.parse.quote(p.get('brand','')),timeout=900))
            for b in r['top']+r['unrated']: b.pop('image',None)
            return p['name'],r
        except Exception as e: err=str(e)
    return p['name'],{'error':err}
with cf.ThreadPoolExecutor(3) as ex: R=dict(ex.map(f,P))
bad=[n for n,r in R.items() if 'error' in r]
import sys
src=sys.argv[1] if len(sys.argv)>1 else 'public/index.html'
out=sys.argv[2] if len(sys.argv)>2 else 'standalone/mezcal_tasting_offline.html'
html=open(src).read()
shim='''<script>
(function(){const D=%s;const K='mezcalPalenques';const rf=window.fetch.bind(window);
const get=()=>{try{const s=localStorage.getItem(K);if(s)return JSON.parse(s)}catch(e){}return D.palenques};
window.fetch=async function(u,o){u=String(u);
 if(u.startsWith('/api/palenques')){if(o&&o.method==='PUT'){try{localStorage.setItem(K,o.body)}catch(e){}window.__p=JSON.parse(o.body);return new Response('{"ok":true}')}return new Response(JSON.stringify(window.__p||get()))}
 if(u.startsWith('/api/mezcals')){const n=new URL(u,'http://x').searchParams.get('name');return new Response(JSON.stringify(D.results[n]||{top:[],unrated:[],candidates:0,source:'offline',matched:[]}))}
 return rf(u,o)}})();
</script>
<style>.bottle img{display:none}</style>'''%json.dumps({'palenques':P,'results':R},ensure_ascii=False).replace('</','<\\/')
html=html.replace('<script>',shim+'\n<script>',1)
assert 'mezcalPalenques' in html
os.makedirs('standalone',exist_ok=True)
open(out,'w').write(html)
print('built',len(html)//1024,'KB; failed:',bad)
