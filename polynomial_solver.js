// polynomial_solver.js
/* eslint-disable no-console */
const fs = require('fs');

// -------------------- Utilities --------------------

function gcdBigInt(a, b) {
  a = a < 0n ? -a : a;
  b = b < 0n ? -b : b;
  while (b !== 0n) {
    const t = b;
    b = a % b;
    a = t;
  }
  return a;
}

function simplifyFrac(num, den) {
  if (den === 0n) throw new Error('Zero denominator in fraction');
  if (den < 0n) { num = -num; den = -den; }
  const g = gcdBigInt(num < 0n ? -num : num, den);
  return [num / g, den / g];
}

function addFrac(n1, d1, n2, d2) {
  // (n1/d1) + (n2/d2) = (n1*d2 + n2*d1) / (d1*d2)
  const num = n1 * d2 + n2 * d1;
  const den = d1 * d2;
  // reduce occasionally to keep numbers manageable
  return simplifyFrac(num, den);
}

function mulFracBigInt(n, d, k) {
  return simplifyFrac(n * k, d);
}


function binomBigInt(n, k) {
  if (k < 0 || k > n) return 0n;
  k = Math.min(k, n - k);
  let num = 1n, den = 1n;
  for (let i = 1n; i <= BigInt(k); i++) {
    num *= BigInt(n) - (i - 1n);
    den *= i;
    const g = gcdBigInt(num, den);
    num /= g; den /= g;
  }
  return num / den; 
}



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



function lagrangeConstantAtZero(points, k) {
  const selected = points.slice(0, k);
  if (selected.length < k) throw new Error(`Need at least ${k} points to interpolate`);

  // Validate duplicate x
  for (let i = 0; i < selected.length; i++) {
    for (let j = i + 1; j < selected.length; j++) {
      if (selected[i][0] === selected[j][0]) {
        throw new Error(`Duplicate x detected at x=${selected[i][0]}`);
      }
    }
  }


  const isConsecutive = selected.every(([x], idx) => x === idx + 1);
  if (isConsecutive) {
    let sum = 0n;
    for (let i = 0; i < k; i++) {
      const sign = (i % 2 === 0) ? 1n : -1n; // (-1)^(i)
      const coeff = binomBigInt(k, i + 1);   // C(k, i+1)
      const yi = selected[i][1];
      sum += sign * coeff * yi;
    }
    return sum; // already integral
  }


  let numSum = 0n;
  let denSum = 1n;

  for (let i = 0; i < selected.length; i++) {
    const [xi, yi] = selected[i];

    let num = 1n;
    let den = 1n;
    for (let j = 0; j < selected.length; j++) {
      if (j === i) continue;
      const [xj] = selected[j];
      num *= BigInt(0 - xj);
      den *= BigInt(xi - xj);
    }
  
    let [tn, td] = simplifyFrac(num * yi, den);
  
    [numSum, denSum] = addFrac(numSum, denSum, tn, td);
  }

  if (numSum % denSum !== 0n) {
    throw new Error(`Non-integer constant term encountered (num=${numSum}, den=${denSum}); input may be inconsistent or insufficient.`);
  }
  return numSum / denSum;
}


function lagrangeEvaluate(points, k, xEval) {
  const selected = points.slice(0, k);
  let numSum = 0n;
  let denSum = 1n;

  for (let i = 0; i < selected.length; i++) {
    const [xi, yi] = selected[i];
    let num = 1n;
    let den = 1n;
    for (let j = 0; j < selected.length; j++) {
      if (j === i) continue;
      const [xj] = selected[j];
      num *= BigInt(xEval - xj);
      den *= BigInt(xi - xj);
    }
    let [tn, td] = simplifyFrac(num * yi, den);
    [numSum, denSum] = addFrac(numSum, denSum, tn, td);
  }
  if (numSum % denSum !== 0n) {
    // Return a string to aid debugging if it’s fractional
    return `${numSum}/${denSum}`;
  }
  return numSum / denSum;
}

// -------------------- Main solving --------------------

function solvePolynomial(jsonData) {
  const n = jsonData.keys.n;
  const k = jsonData.keys.k;

  console.log(`Number of roots provided: ${n}`);
  console.log(`Minimum roots required: ${k}`);

  const points = [];
  for (let i = 1; i <= n; i++) {
    const key = String(i);
    if (jsonData[key]) {
      const base = parseInt(jsonData[key].base, 10);
      const value = jsonData[key].value;
      const x = i;
      const y = baseToDecimal(value, base);
      points.push([x, y]);
      console.log(`Point ${i}: (${x}, ${y.toString()}) - Base ${base}, Value: ${value}`);
    }
  }

  if (points.length < k) {
    throw new Error(`Only ${points.length} usable points present, but k=${k} required`);
  }

  points.sort((a, b) => a[0] - b[0]);
  const constantTerm = lagrangeConstantAtZero(points, k);
  return { constantTerm, points, k };
}

function solveFromFile(filename) {
  try {
    console.log(`\n=== SOLVING ${filename.toUpperCase()} ===`);
    const data = fs.readFileSync(filename, 'utf8');
    const jsonData = JSON.parse(data);
    const { constantTerm, points, k } = solvePolynomial(jsonData);
    console.log(`\nCONSTANT TERM (c) for ${filename}: ${constantTerm.toString()}`);
    console.log('='.repeat(50));

    // Optional: consistency check
    console.log('Consistency check (does each provided point lie on the degree-(k-1) polynomial from the first k points?):');
    for (const [x, y] of points) {
      const px = lagrangeEvaluate(points, k, x);
      const ok = typeof px === 'bigint' ? (px === y) : false;
      console.log(`x=${x}: P(x)=${px.toString()} vs y=${y.toString()} => ${ok ? 'OK' : 'MISMATCH'}`);
    }
    console.log('='.repeat(50));

    return constantTerm;
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



function verifyTestCase1() {
  console.log('\n=== MANUAL VERIFICATION FOR TEST CASE 1 ===');
  console.log('Points after base conversion: (1,4), (2,7), (3,12), (6,39)');
  console.log('Using first 3 points gives P(x) = x^2 + 3, so constant term c = 3');
}

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
  lagrangeConstantAtZero,
  lagrangeEvaluate,
  createTestFiles
};

if (require.main === module) {
  main();
}
