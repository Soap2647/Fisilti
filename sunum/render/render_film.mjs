// Fısıltı filmini kare kare yakalayıp ffmpeg ile MP4'e çevirir.
// Kullanım:
//   node render_film.mjs                 -> ../Fisilti_Film.mp4 (tam film, 1920x1080, 30 fps)
//   node render_film.mjs --test 5,40,120 -> belirtilen saniyelerden test kareleri (PNG)
import { spawn } from 'node:child_process';
import { mkdirSync } from 'node:fs';
import { dirname, resolve } from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';
import puppeteer from 'puppeteer-core';

const here = dirname(fileURLToPath(import.meta.url));
const page_url = pathToFileURL(resolve(here, '..', 'fisilti_film.html')).href + '#render';
const CHROME = 'C:\\Program Files\\Google\\Chrome\\Application\\chrome.exe';
const W = 1920, H = 1080;
const args = process.argv.slice(2);
const testIdx = args.indexOf('--test');
const testSeconds = testIdx >= 0 ? args[testIdx + 1].split(',').map(Number) : null;

const browser = await puppeteer.launch({
  executablePath: CHROME,
  headless: 'new',
  args: ['--ignore-gpu-blocklist', '--enable-webgl', '--use-angle=swiftshader', '--enable-unsafe-swiftshader',
         `--window-size=${W},${H}`, '--hide-scrollbars'],
});
const page = await browser.newPage();
await page.setViewport({ width: W, height: H, deviceScaleFactor: 1 });
page.on('pageerror', e => console.error('[sayfa hatası]', e.message));
await page.goto(page_url, { waitUntil: 'networkidle0', timeout: 120000 });
await page.waitForFunction(() => window.__film && document.fonts.status === 'loaded', { timeout: 60000 });
const { fps, total } = await page.evaluate(() => ({ fps: window.__film.fps, total: window.__film.total }));
const totalFrames = Math.round(total * fps);
const nextPaint = () => page.evaluate(() => new Promise(r => requestAnimationFrame(() => r())));

if (testSeconds) {
  mkdirSync(resolve(here, 'test'), { recursive: true });
  let f = 0;
  for (const sec of testSeconds.sort((a, b) => a - b)) {
    const target = Math.round(sec * fps);
    while (f < target) { await page.evaluate(() => window.__film.step()); f++; }
    await nextPaint();
    const out = resolve(here, 'test', `kare_${String(sec).padStart(3, '0')}s.png`);
    await page.screenshot({ path: out });
    console.log('[ok] test karesi', out);
  }
  await browser.close();
  process.exit(0);
}

const outFile = resolve(here, '..', 'Fisilti_Film.mp4');
const ff = spawn('ffmpeg', ['-y', '-loglevel', 'error', '-f', 'image2pipe', '-framerate', String(fps), '-c:v', 'mjpeg', '-i', '-',
  '-c:v', 'libx264', '-preset', 'slow', '-crf', '18', '-pix_fmt', 'yuv420p', '-movflags', '+faststart', outFile],
  { stdio: ['pipe', 'inherit', 'inherit'] });
const t0 = Date.now();
for (let f = 0; f < totalFrames; f++) {
  await page.evaluate(() => window.__film.step());
  await nextPaint();
  const buf = await page.screenshot({ type: 'jpeg', quality: 95 });
  if (!ff.stdin.write(buf)) await new Promise(r => ff.stdin.once('drain', r));
  if (f % 150 === 0) {
    const el = (Date.now() - t0) / 1000, eta = el / (f + 1) * (totalFrames - f - 1);
    console.log(`[..] kare ${f}/${totalFrames}  geçen ${el.toFixed(0)} sn  kalan ~${(eta / 60).toFixed(1)} dk`);
  }
}
ff.stdin.end();
await new Promise((res, rej) => ff.on('close', c => (c === 0 ? res() : rej(new Error('ffmpeg çıkış kodu ' + c)))));
await browser.close();
console.log(`[ok] ${outFile} · ${totalFrames} kare · ${total} sn · ${((Date.now() - t0) / 60000).toFixed(1)} dk`);
