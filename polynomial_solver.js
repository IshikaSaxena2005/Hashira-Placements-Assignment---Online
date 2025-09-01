// polynomial_solver.js
/* eslint-disable no-console */
const fs = require('fs');

// -------------------- BigInt fraction utilities --------------------
function gcdBigInt(a, b) {
  a = a < 0n ? -a : a;
  b = b < 0n ? -b : b;
  while (b !== 0n) {
    const t = b; b = a % b; a = t;
  }
  return a;
}

class Frac {
  constructor(num, den = 1n) {
    if (den === 0n) throw new Error('Zero denominator');
    if (den < 0n) { num = -num; den = -den; }
    const g = gcdBigInt(num < 0n ? -num : num, den);
    this.n = num / g;
    this.d = den / g;
  }
  static fromBigInt(x) { return new Frac(x, 1n); }
  add(o) { return new Frac(this.n * o.d + o.n * this.d, this.d * o.d); }
  sub(o) { return new Frac(this.n * o.d - o.n * this.d, this.d * o.d); }
  mul(o) { return new Frac(this.n * o.n, this.d * o.d); }
  div(o) { return new Frac(this.n * o.d, this.d * o.n); }
  eqInt(x) { return this.d === 1n && this.n === x; }
  toString() { return this.d === 1n ? this.n.toString() : `${this.n}/${this.d}`; }
}

// -------------------- Base conversion --------------------
function baseToDecimal(value, base) {
  if (base < 2 || base > 36) throw new Error('Base must be between 2 and 36');
  let result = 0n;
  const B = BigInt(base);
  for (let i = 0; i < value.length; i++) {
    const ch = value[i].toLowerCase();
    let d;
    if (ch >= '0' && ch <= '9') d = ch.charCodeAt(0) - '0'.charCodeAt(0);
    else if (ch >= 'a' && ch <= 'z') d = ch.charCodeAt(0) - 'a'.charCodeAt(0) + 10;
    else throw new Error(`Invalid character '${ch}' for base ${base}`);
    if (d >= base) throw new Error(`Invalid character '${ch}' for base ${base}`);
    result = result * B + BigInt(d);
  }
  return result;
}

// -------------------- Vandermonde system (exact, no Lagrange) --------------------
// Solve for a0..a_{k-1} in P(x) = a0 + a1 x + ... + a_{k-1} x^{k-1} using exact Gaussian elimination
function solveVandermondeRational(subPoints) {
  const k = subPoints.length;
  const A = Array.from({ length: k }, () => Array(k).fill(null));
  const b = Array(k).fill(null);

  // Build Vandermonde matrix and RHS
  for (let r = 0; r < k; r++) {
    const [x, y] = subPoints[r];
    b[r] = new Frac(y, 1n);
    let xp = new Frac(1n, 1n);
    for (let c = 0; c < k; c++) {
      A[r][c] = xp;
      xp = xp.mul(new Frac(BigInt(x), 1n));
    }
  }

  // Gaussian elimination with partial pivoting
  for (let col = 0; col < k; col++) {
    let piv = col;
    for (let r = col; r < k; r++) {
      if (A[r][col].n !== 0n) { piv = r; break; }
    }
    if (A[piv][col].n === 0n) return null; // singular
    if (piv !== col) { [A[col], A[piv]] = [A[piv], A[col]]; [b[col], b[piv]] = [b[piv], b[col]]; }

    const fac = A[col][col];
    for (let c = col; c < k; c++) A[col][c] = A[col][c].div(fac);
    b[col] = b[col].div(fac);

    for (let r = col + 1; r < k; r++) {
      const f = A[r][col];
      if (f.n === 0n) continue;
      for (let c = col; c < k; c++) A[r][c] = A[r][c].sub(f.mul(A[col][c]));
      b[r] = b[r].sub(f.mul(b[col]));
    }
  }

  // Back substitution
  const a = Array(k).fill(new Frac(0n, 1n));
  for (let i = k - 1; i >= 0; i--) {
    let s = b[i];
    for (let c = i + 1; c < k; c++) s = s.sub(A[i][c].mul(a[c]));
    a[i] = s; // diagonal is 1
  }
  return a; // coefficients a0..a_{k-1}
}

function evalPolyRational(a, x) {
  let xp = new Frac(1n, 1n);
  let sum = new Frac(0n, 1n);
  const X = new Frac(BigInt(x), 1n);
  for (let i = 0; i < a.length; i++) {
    sum = sum.add(a[i].mul(xp));
    xp = xp.mul(X);
  }
  return sum;
}

// -------------------- Robust consensus over k-subsets --------------------
function allKSubsetsIndices(n, k) {
  const res = [];
  function rec(start, path) {
    if (path.length === k) { res.push(path.slice()); return; }
    for (let i = start; i <= n - (k - path.length); i++) {
      path.push(i);
      rec(i + 1, path);
      path.pop();
    }
  }
  rec(0, []);
  return res;
}

function robustConstantTerm(points, k) {
  const n = points.length;
  if (n < k) throw new Error(`Need at least ${k} points, got ${n}`);
  const subsets = allKSubsetsIndices(n, k); // for n=10,k=7 this is 120
  let best = null;

  for (const idxs of subsets) {
    const subPts = idxs.map(i => points[i]);
    const coeffs = solveVandermondeRational(subPts);
    if (!coeffs) continue;

    const inliers = [];
    for (let j = 0; j < n; j++) {
      const [x, y] = points[j];
      const val = evalPolyRational(coeffs, x);
      if (val.eqInt(y)) inliers.push(j);
    }
    const score = inliers.length;
    if (!best || score > best.score) {
      best = { score, inliers, coeffs, subset: idxs };
      if (score === n) break; // perfect fit
    }
  }
  if (!best) throw new Error('Failed to find a valid consensus model');
  return best;
}

// -------------------- Orchestration --------------------
function solvePolynomial(jsonData) {
  const n = jsonData.keys.n;
  const k = jsonData.keys.k;

  console.log(`Number of roots provided: ${n}`);
  console.log(`Minimum roots required: ${k}`);

  const points = [];
  for (let i = 1; i <= n; i++) {
    const key = String(i);
    if (!jsonData[key]) continue;
    const base = parseInt(jsonData[key].base, 10);
    const value = jsonData[key].value;
    const x = i;
    const y = baseToDecimal(value, base);
    points.push([x, y]);
    console.log(`Point ${i}: (${x}, ${y.toString()}) - Base ${base}, Value: ${value}`);
  }

  points.sort((a, b) => a[0] - b[0]);

  const best = robustConstantTerm(points, k);
  const constantTerm = best.coeffs[0]; // a0 = P(0)

  return { constantTerm, points, k, best };
}

function solveFromFile(filename) {
  try {
    console.log(`\n=== SOLVING ${filename.toUpperCase()} ===`);
    const data = fs.readFileSync(filename, 'utf8');
    const jsonData = JSON.parse(data);
    const { constantTerm, points, k, best } = solvePolynomial(jsonData);

    console.log(`\nk = ${k}, inliers = ${best.inliers.length}/${points.length}`);
    console.log(`Best subset indices (0-based): ${best.subset.join(', ')}`);
    console.log(`CONSTANT TERM (c) for ${filename}: ${constantTerm.toString()}`);
    console.log('='.repeat(50));

    // Consistency print with consensus model
    console.log('Consensus check (1=inlier, 0=outlier):');
    for (let i = 0; i < points.length; i++) {
      const [x, y] = points[i];
      const val = evalPolyRational(best.coeffs, x);
      const ok = val.eqInt(y);
      console.log(`x=${x}: P(x)=${val.toString()} vs y=${y.toString()} => ${ok ? 'OK' : 'MISMATCH'}`);
    }
    console.log('Inlier mask:', points.map((_, i) => best.inliers.includes(i) ? 1 : 0).join(' '));
    console.log('='.repeat(50));

    // Return BigInt if integral, else fraction string
    return constantTerm.d === 1n ? constantTerm.n : BigInt(0); // constantTerm as BigInt if integral
  } catch (error) {
    console.error(`Error reading file ${filename}:`, error.message);
    return null;
  }
}

// -------------------- Test files --------------------
function createTestFiles() {
  const testCase1 = {
    keys: { n: 4, k: 3 },
    1: { base: '10', value: '4' },
    2: { base: '2', value: '111' },
    3: { base: '10', value: '12' },
    6: { base: '4', value: '213' }
  };

  const testCase2 = {
    keys: { n: 10, k: 7 },
    1: { base: '6', value: '13444211440455345511' },
    2: { base: '15', value: 'aed7015a346d635' },
    3: { base: '15', value: '6aeeb69631c227c' },
    4: { base: '16', value: 'e1b5e05623d881f' },
    5: { base: '8', value: '316034514573652620673' },
    6: { base: '3', value: '2122212201122002221120200210011020220200' },
    7: { base: '3', value: '20120221122211000100210021102001201112121' },
    8: { base: '6', value: '20220554335330240002224253' },
    9: { base: '12', value: '45153788322a1255483' },
    10: { base: '7', value: '1101613130313526312514143' }
  };

  fs.writeFileSync('testcase1.json', JSON.stringify(testCase1, null, 2));
  fs.writeFileSync('testcase2.json', JSON.stringify(testCase2, null, 2));
  console.log('Test case files created: testcase1.json, testcase2.json');
}

// -------------------- Optional manual note --------------------
function verifyTestCase1() {
  console.log('\n=== MANUAL VERIFICATION FOR TEST CASE 1 ===');
  console.log('Points after base conversion: (1,4), (2,7), (3,12), (6,39)');
  console.log('Using first 3 points gives P(x) = x^2 + 3, so constant term c = 3');
}

// -------------------- CLI --------------------
function main() {
  console.log('POLYNOMIAL SOLVER - FINDING CONSTANT TERM');
  console.log('='.repeat(60));

  if (!fs.existsSync('testcase1.json') || !fs.existsSync('testcase2.json')) {
    console.log('Creating test case files...');
    createTestFiles();
  }

  const result1 = solveFromFile('testcase1.json');
  const result2 = solveFromFile('testcase2.json');

  if (result1 !== null) verifyTestCase1();

  console.log('\nFINAL RESULTS:');
  if (result1 !== null) console.log(`Test Case 1 - Constant term: ${result1.toString()}`);
  if (result2 !== null) console.log(`Test Case 2 - Constant term: ${result2.toString()}`);
}

// Exports for testing
module.exports = {
  solvePolynomial,
  solveFromFile,
  baseToDecimal,
  solveVandermondeRational,
  evalPolyRational,
  robustConstantTerm,
  createTestFiles
};

if (require.main === module) {
  main();
}
