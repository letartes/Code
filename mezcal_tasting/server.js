// Mexico in a Bottle tasting planner. Zero dependencies, Node 18+.
// Run: node server.js  then open http://localhost:3055
const http = require('http');
const fs = require('fs');
const path = require('path');

const PORT = process.env.PORT || 3055;
const BASE = 'https://www.mezcalreviews.com';
const DATA = path.join(__dirname, 'data');
const PALENQUES = path.join(DATA, 'palenques.json');
const CACHE = path.join(DATA, 'cache.json');
const UA = 'Mozilla/5.0 (Macintosh) AppleWebKit/537.36 Chrome/124 Safari/537.36';
const CACHE_TTL = 7 * 24 * 3600 * 1000;

const readJSON = (f, d) => { try { return JSON.parse(fs.readFileSync(f, 'utf8')); } catch { return d; } };
const writeJSON = (f, o) => fs.writeFileSync(f, JSON.stringify(o, null, 2));
const slugify = s => s.toLowerCase().normalize('NFD').replace(/[̀-ͯ]/g, '')
  .replace(/[^a-z0-9]+/g, '-').replace(/^-|-$/g, '');

async function get(url) {
  for (let i = 0; i < 3; i++) {
    try {
      const r = await fetch(url, { headers: { 'User-Agent': UA }, redirect: 'follow' });
      if (r.ok) return r.text();
      if (r.status === 404) return null;
    } catch {}
    await new Promise(res => setTimeout(res, 600 * (i + 1)));
  }
  return null;
}
const text = html => html.replace(/<script[\s\S]*?<\/script>|<style[\s\S]*?<\/style>/g, ' ')
  .replace(/<[^>]+>/g, ' ').replace(/&nbsp;|&#160;/g, ' ').replace(/\s+/g, ' ');
const decode = s => s.replace(/&#8211;|&ndash;/g, '-').replace(/&#8217;|&rsquo;/g, "'")
  .replace(/&amp;/g, '&').replace(/&#038;/g, '&').replace(/&#(\d+);/g, (_, n) => String.fromCharCode(n));

function mezcalLinks(html) {
  const out = new Map();
  // Only listing cards; sidebar widgets also link to unrelated "top rated" bottles.
  for (const m of html.matchAll(/<a href="(https?:\/\/(?:www\.)?mezcalreviews\.com\/mezcal\/([a-z0-9-]+)\/?)"\s+class="mr__review-box__inner/gi))
    out.set(m[2], BASE + '/mezcal/' + m[2] + '/');
  return [...out.entries()].map(([slug, url]) => ({ slug, url }));
}

function parseBottle(html, url) {
  const t = text(html);
  let rating = null, reviews = 0;
  let m = t.match(/([\d.]+)\s*out of\s*5\s*based on\s*(\d+)\s*reviews?/i);
  if (m) { rating = parseFloat(m[1]); reviews = parseInt(m[2], 10); }
  else if ((m = html.match(/"ratingValue"\s*:\s*"?([\d.]+)"?[\s\S]{0,200}?"(?:ratingCount|reviewCount)"\s*:\s*"?(\d+)/))) {
    rating = parseFloat(m[1]); reviews = parseInt(m[2], 10);
  }
  const h1 = html.match(/<h1[^>]*>([\s\S]*?)<\/h1>/i);
  const title = h1 ? decode(text(h1[1]).trim()) : url.split('/').filter(Boolean).pop();
  const desc = html.match(/<meta property="og:description" content="([^"]*)"/i);
  const img = html.match(/<meta property="og:image" content="([^"]*)"/i);
  return { title, url, rating, reviews, blurb: desc ? decode(desc[1]).slice(0, 220) : '', image: img ? img[1] : '' };
}

// Weighted score so a lone 5.0 does not beat a 4.5 with many reviews.
const score = b => b.rating == null ? -1 : (b.reviews * b.rating + 3 * 3.5) / (b.reviews + 3);

// ---- Fuzzy brand matching -------------------------------------------------
const STOP = new Set(['mezcal','mezcales','spirits','spirit','tequila','raicilla','bacanora','sotol','gin','licor',
  'liqueurs','liqueur','bitters','aguardiente','vodka','rum','whiskey','brewhouse','and','de','del','la','las','los','el','y']);
const NUM = { uno: '1', one: '1', dos: '2', two: '2', tres: '3', three: '3', cuatro: '4', four: '4', cinco: '5', five: '5',
  seis: '6', six: '6', siete: '7', seven: '7', ocho: '8', eight: '8', nueve: '9', nine: '9', diez: '10', ten: '10' };
const norm = s => slugify(s).split('-').filter(w => w && !STOP.has(w)).map(w => NUM[w] || w);
const compact = s => norm(s).join('');
function lev(a, b) {
  const d = Array.from({ length: a.length + 1 }, (_, i) => [i, ...Array(b.length).fill(0)]);
  for (let j = 1; j <= b.length; j++) d[0][j] = j;
  for (let i = 1; i <= a.length; i++) for (let j = 1; j <= b.length; j++)
    d[i][j] = Math.min(d[i - 1][j] + 1, d[i][j - 1] + 1, d[i - 1][j - 1] + (a[i - 1] === b[j - 1] ? 0 : 1));
  return d[a.length][b.length];
}
function sim(x, y) {
  const a = compact(x), b = compact(y);
  if (!a || !b) return 0;
  if (a === b) return 1;
  let s = 1 - lev(a, b) / Math.max(a.length, b.length);
  const [sh, lg] = a.length < b.length ? [a, b] : [b, a];
  if (sh.length >= 5 && lg.includes(sh)) s = Math.max(s, 0.8 + 0.2 * sh.length / lg.length);
  return s;
}

async function brandList() {
  const cache = readJSON(CACHE, {});
  if (cache.__brands && Date.now() - cache.__brands.at < 30 * 24 * 3600 * 1000) return cache.__brands.list;
  // Each brand page's dropdown omits its own brand, so union the lists from two different pages.
  const map = new Map();
  for (const seed of ['real-minero', 'alipus']) {
    const html = await get(`${BASE}/filter-by/brand/${seed}/`);
    if (!html) continue;
    for (const m of html.matchAll(/<option value="https?:\/\/[^"]*\/filter-by\/brand\/([^"\/]+)\/">([^<]*)<\/option>/g)) {
      const lab = decode(m[2]); const c = lab.match(/\((\d+)\)\s*$/);
      map.set(m[1], { slug: m[1], label: lab.replace(/\s*\(\d+\)\s*$/, ''), count: c ? +c[1] : 0 });
    }
  }
  const list = [...map.values()];
  if (list.length) { cache.__brands = { at: Date.now(), list }; writeJSON(CACHE, cache); }
  return list;
}

async function brandPages(slug, max = 10) {
  const links = [];
  for (let p = 1; p <= max; p++) {
    const html = await get(`${BASE}/filter-by/brand/${slug}/${p > 1 ? 'page/' + p + '/' : ''}`);
    if (!html) break;
    const found = mezcalLinks(html).filter(l => !links.some(x => x.slug === l.slug));
    if (!found.length) break;
    links.push(...found);
  }
  return links;
}

async function findCandidates(name, brandSlug) {
  let links = [], matched = [], source = '';
  if (brandSlug) {
    matched = [brandSlug]; links = await brandPages(brandSlug); source = 'brand archive';
  } else {
    const brands = (await brandList()).map(b => ({ ...b, s: sim(name, b.label) })).sort((a, b) => b.s - a.s);
    const best = brands[0];
    if (best && best.s >= 0.8) {
      const picks = brands.filter(b => b.s >= 0.8 && b.s >= best.s - 0.05).slice(0, 2);
      for (const b of picks) { matched.push(b.label); links.push(...await brandPages(b.slug)); }
      source = 'brand archive';
    }
  }
  if (!links.length) {
    // Fallback: site search, keep bottles whose title shares a distinctive word with the name.
    source = 'site search';
    const words = norm(name).filter(w => w.length > 2);
    const html = await get(`${BASE}/?s=${encodeURIComponent(words.join(' ') || name)}`);
    if (html) links = mezcalLinks(html).filter(l => {
      const sw = l.slug.split('-');
      return words.some(w => sw.some(x => x === w || (w.length > 4 && sim(w, x) >= 0.82)));
    });
  }
  const seen = new Set();
  return { source, matched, links: links.filter(l => !seen.has(l.slug) && seen.add(l.slug)).slice(0, 120) };
}

async function topMezcals(name, brandSlug, force) {
  const key = slugify(name) + '|' + (brandSlug || '');
  const cache = readJSON(CACHE, {});
  if (!force && cache[key] && Date.now() - cache[key].at < CACHE_TTL) return cache[key];
  const { source, matched, links } = await findCandidates(name, brandSlug);
  const bottles = [];
  for (let i = 0; i < links.length; i += 6) {
    const batch = await Promise.all(links.slice(i, i + 6).map(async l => {
      try { const h = await get(l.url); return h ? parseBottle(h, l.url) : null; } catch { return null; }
    }));
    bottles.push(...batch.filter(Boolean));
  }
  const ranked = bottles.filter(b => b.rating != null).sort((a, b) => score(b) - score(a));
  const unrated = bottles.filter(b => b.rating == null);
  const result = {
    at: Date.now(), source, matched, candidates: links.length,
    top: ranked.map(b => ({ ...b, score: +score(b).toFixed(2) })),
    unrated,
  };
  if (links.length) { cache[key] = result; writeJSON(CACHE, cache); }
  return result;
}

const types = { '.html': 'text/html', '.js': 'text/javascript', '.css': 'text/css' };
http.createServer(async (req, res) => {
  const u = new URL(req.url, 'http://x');
  const send = (code, body, type = 'application/json') => { res.writeHead(code, { 'Content-Type': type }); res.end(typeof body === 'string' || Buffer.isBuffer(body) ? body : JSON.stringify(body)); };
  try {
    if (u.pathname === '/api/palenques' && req.method === 'GET') return send(200, readJSON(PALENQUES, []));
    if (u.pathname === '/api/palenques' && req.method === 'PUT') {
      let b = ''; for await (const c of req) b += c;
      writeJSON(PALENQUES, JSON.parse(b)); return send(200, { ok: true });
    }
    if (u.pathname === '/api/mezcals') {
      const name = u.searchParams.get('name');
      if (!name) return send(400, { error: 'name required' });
      return send(200, await topMezcals(name, u.searchParams.get('brand') || '', u.searchParams.get('force') === '1'));
    }
    const file = path.join(__dirname, 'public', u.pathname === '/' ? 'index.html' : u.pathname);
    if (file.startsWith(path.join(__dirname, 'public')) && fs.existsSync(file))
      return send(200, fs.readFileSync(file), types[path.extname(file)] || 'application/octet-stream');
    send(404, { error: 'not found' });
  } catch (e) { send(500, { error: String(e.message || e) }); }
}).listen(PORT, () => console.log(`Mezcal tasting planner: http://localhost:${PORT}`));
