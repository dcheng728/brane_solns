import { createServer } from 'node:http';
import { readFileSync, mkdirSync, existsSync, readdirSync, statSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { spawnSync } from 'node:child_process';
import { chromium } from 'playwright';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const repoRoot = path.resolve(__dirname, '..');
const docsDir = path.join(repoRoot, 'docs');
const outDir = path.join(repoRoot, 'build', 'docs-pdf');
const htmlDir = path.join(repoRoot, 'build', 'docs-html');
const katexHeaderPath = path.join(__dirname, 'katex-header.html');

function walk(dir) {
  const entries = readdirSync(dir, { withFileTypes: true });
  const files = [];
  for (const entry of entries) {
    const full = path.join(dir, entry.name);
    if (entry.isDirectory()) {
      if (entry.name === 'images' || entry.name.startsWith('.')) continue;
      files.push(...walk(full));
    } else if (entry.isFile()) {
      if (entry.name.toLowerCase().endsWith('.md')) files.push(full);
    }
  }
  return files;
}

function ensureDir(p) {
  mkdirSync(p, { recursive: true });
}

function runPandoc({ inputMd, outputHtml }) {
  const header = existsSync(katexHeaderPath) ? katexHeaderPath : null;
  const args = [
    '-f', 'gfm',
    '-t', 'html5',
    '-s',
    '--metadata', `title=${path.basename(inputMd, '.md')}`,
    '--resource-path', `${docsDir}${path.delimiter}${repoRoot}`,
  ];

  if (header) {
    args.push('--include-in-header', header);
  }

  args.push('-o', outputHtml, inputMd);

  const result = spawnSync('pandoc', args, {
    cwd: repoRoot,
    encoding: 'utf8',
    shell: false,
  });

  if (result.error) throw result.error;
  if (result.status !== 0) {
    throw new Error(`pandoc failed for ${inputMd}\n${result.stderr || result.stdout || ''}`);
  }
}

function startStaticServer(rootDir) {
  const server = createServer((req, res) => {
    try {
      const urlPath = decodeURIComponent((req.url || '/').split('?')[0]);
      const safePath = urlPath === '/' ? '/index.html' : urlPath;
      const absPath = path.resolve(rootDir, '.' + safePath);

      if (!absPath.startsWith(rootDir)) {
        res.statusCode = 403;
        res.end('Forbidden');
        return;
      }

      if (!existsSync(absPath) || !statSync(absPath).isFile()) {
        res.statusCode = 404;
        res.end('Not Found');
        return;
      }

      const buf = readFileSync(absPath);
      const ext = path.extname(absPath).toLowerCase();
      const contentType =
        ext === '.html' ? 'text/html; charset=utf-8' :
        ext === '.css' ? 'text/css; charset=utf-8' :
        ext === '.js' ? 'text/javascript; charset=utf-8' :
        ext === '.svg' ? 'image/svg+xml' :
        ext === '.png' ? 'image/png' :
        ext === '.jpg' || ext === '.jpeg' ? 'image/jpeg' :
        ext === '.gif' ? 'image/gif' :
        'application/octet-stream';

      res.setHeader('Content-Type', contentType);
      res.end(buf);
    } catch (e) {
      res.statusCode = 500;
      res.end(String(e));
    }
  });

  return new Promise((resolve) => {
    server.listen(0, '127.0.0.1', () => {
      const addr = server.address();
      resolve({ server, port: addr.port });
    });
  });
}

async function main() {
  if (!existsSync(docsDir)) {
    throw new Error(`docs directory not found at ${docsDir}`);
  }

  ensureDir(outDir);
  ensureDir(htmlDir);

  const mdFiles = walk(docsDir);
  if (mdFiles.length === 0) {
    console.log('No markdown files found under docs/.');
    return;
  }

  // Generate HTML first (Pandoc), then render via Chromium (Playwright) so KaTeX JS runs.
  const htmlFiles = [];
  for (const mdFile of mdFiles) {
    const rel = path.relative(docsDir, mdFile);
    const baseNoExt = rel.replace(/\.md$/i, '');
    const outHtml = path.join(htmlDir, baseNoExt + '.html');

    ensureDir(path.dirname(outHtml));
    runPandoc({ inputMd: mdFile, outputHtml: outHtml });
    htmlFiles.push({ mdFile, rel, outHtml, baseNoExt });
  }

  const { server, port } = await startStaticServer(repoRoot);

  const browser = await chromium.launch();
  try {
    const page = await browser.newPage();

    for (const { outHtml, baseNoExt, mdFile } of htmlFiles) {
      const relHtmlFromRoot = path.relative(repoRoot, outHtml).split(path.sep).join('/');
      const url = `http://127.0.0.1:${port}/${relHtmlFromRoot}`;
      const outPdf = path.join(outDir, baseNoExt + '.pdf');
      ensureDir(path.dirname(outPdf));

      await page.goto(url, { waitUntil: 'networkidle' });
      await page.waitForFunction(() => window.__katexRendered === true, null, { timeout: 30000 }).catch(() => undefined);

      await page.pdf({
        path: outPdf,
        format: 'A4',
        printBackground: true,
        margin: { top: '20mm', right: '18mm', bottom: '20mm', left: '18mm' },
      });

      console.log(`PDF: ${outPdf}  (from ${mdFile})`);
    }

    await page.close();
  } finally {
    await browser.close();
    server.close();
  }
}

main().catch((err) => {
  console.error(err);
  process.exit(1);
});
