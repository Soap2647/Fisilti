// Fısıltı kapanış filmini render eder: görüntü (kare kare) + ses olayları listesi.
// Kullanım:
//   node render_kapanis.mjs --test 5,12,19,25,38,48,58,68  -> test kareleri (PNG)
//   node render_kapanis.mjs                               -> kapanis_goruntu.mp4 + kapanis_cues.json
import { spawn } from 'node:child_process';
import { mkdirSync, writeFileSync } from 'node:fs';
import { dirname, resolve } from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';
import puppeteer from 'puppeteer-core';

const here = dirname(fileURLToPath(import.meta.url));
const url = pathToFileURL(resolve(here, '..', 'fisilti_kapanis.html')).href + '#render';
const CHROME = 'C:\\Program Files\\Google\\Chrome\\Application\\chrome.exe';
const W = 1920, H = 1080;
const args = process.argv.slice(2);
const ti = args.indexOf('--test');
const tests = ti >= 0 ? args[ti + 1].split(',').map(Number) : null;

const browser = await puppeteer.launch({
  executablePath: CHROME, headless: 'new',
  args: ['--ignore-gpu-blocklist', '--enable-webgl', '--use-angle=swiftshader', '--enable-unsafe-swiftshader', `--window-size=${W},${H}`, '--hide-scrollbars'],
});
const page = await browser.newPage();
await page.setViewport({ width: W, height: H, deviceScaleFactor: 1 });
page.on('pageerror', e => console.error('[sayfa hatası]', e.message));
page.on('console', m => { if (m.type() === 'error') console.error('[konsol]', m.text()); });
await page.goto(url, { waitUntil: 'networkidle0', timeout: 120000 });
await page.waitForFunction(() => window.__film && window.__film.ready, { timeout: 60000 });
const { fps, total, cues } = await page.evaluate(() => ({ fps: window.__film.fps, total: window.__film.total, cues: window.__film.cues }));
const paint = () => page.evaluate(() => new Promise(r => requestAnimationFrame(() => r())));

if (tests) {
  mkdirSync(resolve(here, 'test'), { recursive: true });
  for (const s of tests) {
    await page.evaluate(t => window.__film.seek(t), s); await paint();
    const out = resolve(here, 'test', `kapanis_${String(s).padStart(3, '0')}s.png`);
    await page.screenshot({ path: out }); console.log('[ok]', out);
  }
  await browser.close(); process.exit(0);
}

writeFileSync(resolve(here, 'kapanis_cues.json'), JSON.stringify({ total, cues }, null, 1));
const video = resolve(here, 'kapanis_goruntu.mp4');
const ff = spawn('ffmpeg', ['-y', '-loglevel', 'error', '-f', 'image2pipe', '-framerate', String(fps), '-c:v', 'mjpeg', '-i', '-',
  '-c:v', 'libx264', '-preset', 'slow', '-crf', '17', '-pix_fmt', 'yuv420p', '-movflags', '+faststart', video], { stdio: ['pipe', 'inherit', 'inherit'] });
const frames = Math.round(total * fps), t0 = Date.now();
for (let f = 0; f < frames; f++) {
  await page.evaluate(t => window.__film.seek(t), f / fps); await paint();
  const buf = await page.screenshot({ type: 'jpeg', quality: 95 });
  if (!ff.stdin.write(buf)) await new Promise(r => ff.stdin.once('drain', r));
  if (f % 150 === 0) {
    const el = (Date.now() - t0) / 1000;
    console.log(`[..] kare ${f}/${frames}  geçen ${el.toFixed(0)} sn  kalan ~${(el / (f + 1) * (frames - f - 1) / 60).toFixed(1)} dk`);
  }
}
ff.stdin.end();
await new Promise((res, rej) => ff.on('close', c => (c === 0 ? res() : rej(new Error('ffmpeg ' + c)))));
await browser.close();
console.log(`[ok] ${video} · ${frames} kare · ${((Date.now() - t0) / 60000).toFixed(1)} dk`);
