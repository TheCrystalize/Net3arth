let promises = [];
let frame = 0;
let frames = 1;
let render = false;

function sleep(ms) {
  return new Promise(resolve => setTimeout(resolve, ms));
}

async function resolvePromises() {
  while(promises.length > 0) {
    await sleep(1000);
  }
}

/*buffers*/
//default buffer
function buffer(oldBuffer, newBuffer) {
  return {
    red: oldBuffer.red + newBuffer.red * newBuffer.alpha,
    green: oldBuffer.green + newBuffer.green * newBuffer.alpha,
    blue: oldBuffer.blue + newBuffer.blue * newBuffer.alpha,
    z: 0
  }
}

function defaultBuffer(oldBuffer, newBuffer) {
  return {
    red: oldBuffer.red + newBuffer.red * newBuffer.alpha,
    green: oldBuffer.green + newBuffer.green * newBuffer.alpha,
    blue: oldBuffer.blue + newBuffer.blue * newBuffer.alpha,
    z: 0
  }
}

function firstBuffer(oldBuffer, newBuffer) {
  if(oldBuffer.z !== 0) {
    return oldBuffer;
  }
  return {
    red: newBuffer.red * newBuffer.alpha,
    green: newBuffer.green * newBuffer.alpha,
    blue: newBuffer.blue * newBuffer.alpha,
    z: 1
  }
}

function lastBuffer(oldBuffer, newBuffer) {
  return {
    red: newBuffer.red * newBuffer.alpha,
    green: newBuffer.green * newBuffer.alpha,
    blue: newBuffer.blue * newBuffer.alpha,
    z: 1
  }
}

function averageBuffer(oldBuffer, newBuffer) {
  if(newBuffer.alpha === 0) {
    return oldBuffer;
  }

  let total = oldBuffer.z + newBuffer.alpha;
  return {
    red: (oldBuffer.red * oldBuffer.z + newBuffer.red * newBuffer.alpha) / total,
    green: (oldBuffer.green * oldBuffer.z + newBuffer.green * newBuffer.alpha) / total,
    blue: (oldBuffer.blue * oldBuffer.z + newBuffer.blue * newBuffer.alpha) / total,
    z: total,
  }
}

function zBuffer(oldBuffer, newBuffer) {
  if(oldBuffer.z !== 0 && oldBuffer.z < newBuffer.z) {
    return oldBuffer;
  }
  return {
    red: newBuffer.red * newBuffer.alpha,
    green: newBuffer.green * newBuffer.alpha,
    blue: newBuffer.blue * newBuffer.alpha,
    z: newBuffer.z
  };
}

function shaderPass(transform) {
  let newBuffer;
  let newZbuffer;
  return z => {
    if(z.re === z.width - 1 && z.im === z.height - 1) {
      newBuffer = [
        new Float32Array(z.width * z.height),
        new Float32Array(z.width * z.height),
        new Float32Array(z.width * z.height),
      ];
      newZBuffer = new Float64Array(z.width * z.height);
      for(let i = z.zBuffer.length - 1; i >= 0; i--) {
        let ans = transform({
          ...z,
          re: i % z.width,
          im: (i / z.width) >> 0,
          z: z.zBuffer[i],
          red: z.mainBuffer[0][i],
          green: z.mainBuffer[1][i],
          blue: z.mainBuffer[2][i],
        });
        newBuffer[0][i] = ans.red;
        newBuffer[1][i] = ans.green;
        newBuffer[2][i] = ans.blue;
        newZBuffer[i] = ans.z;
      }
    }
    return {
      ...z,
      zBuffer: newZBuffer,
      mainBuffer: newBuffer,
      red: newBuffer[0][z.re + z.im * z.width],
      green: newBuffer[1][z.re + z.im * z.width],
      blue: newBuffer[2][z.re + z.im * z.width],
      z: newZBuffer[z.re + z.im * z.width],
    };
  }
}

function normalizeColors() {
  let brightness;
  return z => {
    if(z.re === z.width - 1 && z.im === z.height - 1) {
      brightness = 0;
      for(let i = z.zBuffer.length - 1; i >= 0; i--) {
        brightness = Math.max(
          brightness,
          z.mainBuffer[0][i],
          z.mainBuffer[1][i],
          z.mainBuffer[2][i]);
      }
    }
    return {
      ...z,
      red: z.red / brightness,
      green: z.green / brightness,
      blue: z.blue / brightness,
    };
  }
}

/* helper functions */
if(!Math.hypot) {
  Math.hypot = function() {
    var y = 0,
      i = arguments.length;
    while(i--) y += arguments[i] * arguments[i];
    return Math.sqrt(y);
  };
}

const DEGREE = Math.PI / 180;

function hash(str) {
  var hash = 0,
    i, chr;
  if (str.length === 0) return hash;
  for (i = 0; i < str.length; i++) {
    chr = str.charCodeAt(i);
    hash = ((hash << 5) - hash) + chr;
    hash |= 0; // Convert to 32bit integer
  }
  return hash;
}

function drawPoint() {}

function draw() {
  return z => {
    drawPoint(z);
    return z;
  };
}

function C(real, imaginary) {
  return {
    re: real,
    im: imaginary
  };
}

function conj(z) {
  return {
    re: z.re,
    im: -z.im
  }
}

function div(z, c) {
  const s = 1 / (c.re * c.re + c.im * c.im);
  return {
    re: (z.re * c.re + z.im * c.im) * s,
    im: (z.im * c.re - z.re * c.im) * s
  }
}

function divScalar(z, n) {
  return {
    re: z.re / n,
    im: z.im / n
  }
}

function add(z, c) {
  return {
    re: z.re + c.re,
    im: z.im + c.im
  }
}

function neg(z) {
  return {
    re: -z.re,
    im: -z.im
  }
}

function sub(z, c) {
  return {
    re: z.re - c.re,
    im: z.im - c.im
  }
}

function addScalar(z, s) {
  return {
    re: z.re + s,
    im: z.im
  }
}

function mult(z, c) {
  return {
    re: z.re * c.re - z.im * c.im,
    im: z.re * c.im + z.im * c.re
  }
}

function multScalar(z, s) {
  return {
    re: z.re * s,
    im: z.im * s
  }
}

function sqrt(z) {
  const s = Math.hypot(z.re, z.im),
    sgn = z.im < 0 ? -1 : 1;
  return multScalar({
    re: Math.sqrt(s + z.re),
    im: sgn * Math.sqrt(s - z.re)
  }, 0.5 * Math.SQRT2);
}

function log(z) {
  return {
    re: 0.5 * Math.log(z.re * z.re + z.im * z.im),
    im: Math.atan2(z.im, z.re)
  }
}

function exp(z) {
  const e = Math.exp(z.re);
  return {
    re: Math.cos(z.im) * e,
    im: Math.sin(z.im) * e
  }
}

function dot(z, c) {
  return (z.re * c.re + z.im * c.im)
}

function sin(z) {
  const im = C(0, 1);
  return add(multScalar(mult(im, exp(mult(neg(im), z))), 0.5), neg(multScalar(mult(im, exp(mult(im, z))), 0.5)))
}

function cos(z) {
  const im = C(0, 1);
  return add(multScalar(exp(mult(neg(im), z)), 0.5), multScalar(exp(mult(im, z)), 0.5))
}

function tan(z) {
  return div(cos(z), sin(z))
}

function sinh(z) {
  return add(divScalar(exp(z), 2), neg(divScalar(exp(neg(z)), 2)))
}

function cosh(z) {
  return add(divScalar(exp(neg(z)), 2), divScalar(exp(z), 2))
}

function tanh(z) {
  return div(sinh(z), cosh(z))
}

function gaussRnd() {
  return (Math.random() + Math.random()) * 2.0 - 2.0;
}

function sign(s) {
  return s < 0 ? -1 : s > 0 ? 1 : 0;
}

function intpow(x, p) {
  let ipret = 1;
  let mulf = x;
  let i;
  let power = Math.floor(p * sign(p));

  if(p < 0) {
    mulf = 1.0 / 1e-8 + Math.abs(mulf);
  }

  for(i = 1; i < power; i++) {
    ipret *= mulf;
  }

  return ipret;
}

function nonz(nz) {
  if(Math.abs(nz) <= Number.EPSILON) {
    return Number.EPSILON * sign(nz);
  }
  return nz;
}

function sacot(angle) {
  if(Math.abs(angle) <= Number.EPSILON) {
    return 0;
  }
  return Math.atan2(1.0, angle);
}

function scot(angle) {
  return Math.cos(angle) / nonz(Math.sin(angle));
}

function stan(angle) {
  return Math.sin(angle) / nonz(Math.cos(angle));
}

function ssqr(s) {
  return s * s;
}

function ssqrt(ssq) {
  return sign(ssq) * Math.sqrt(Math.abs(ssq));
}

function ssqrt2(ssq) {
  return Math.sqrt(Math.abs(ssq));
}

function carlsonRF(x, y, z) {
  let result = 0;
  let a = 0;
  let lambda = 0;
  let dx = 0;
  let dy = 0;
  let dz = 0;
  let minError = 1e-5;
  let mIt = 25;
  let itC = 0;

  do {
    lambda = Math.sqrt(x * y) + Math.sqrt(y * z) + Math.sqrt(z * x);
    x = (x + lambda) * 0.25;
    y = (y + lambda) * 0.25;
    z = (z + lambda) * 0.25;

    a = (x + y + z) / 3;

    dx = 1 - x / a;
    dy = 1 - y / a;
    dz = 1 - z / a;
    itC++;
    if(itC >= mIt) {
      break;
    }
  } while(Math.max(Math.max(Math.abs(dx), Math.abs(dy)), Math.abs(dz)) > minError);

  let e2 = dx * dy + dy * dz + dz * dx;
  let e3 = dy * dx * dz;

  result = 1 - 0.1 * e2 + (1 / 14) * e3 + (1 / 24) * intpow(e2, 2) - (3 / 44) * e2 * e3 - (5 / 208) * intpow(e2, 3) + (3 / 104) * intpow(e3, 2) + (0.0625) * intpow(e2, 2) * e3;
  result *= (1 / Math.sqrt(a));
  return result;
}

function carlsonRD(x, y, z) {
  let alamb, ave, delx, dely, delz, ea, eb, ec, ed, ee, fac, sqrtx, sqrty, sqrtz, sum, xt, yt, zt;
  let c1 = 3 / 14,
    c2 = 1 / 6,
    c3 = 9 / 22,
    c4 = 3 / 26,
    c5 = 0.25 * c3,
    c6 = 1.5 * c4;
  let minError = 1e-5;
  let mIt = 25;
  let itC = 0;
  xt = x;
  yt = y;
  zt = z;
  sum = 0;
  fac = 1;
  do {
    sqrtx = xt ** 0.5;
    sqrty = yt ** 0.5;
    sqrtz = zt ** 0.5;
    alamb = sqrtz * (sqrty + sqrtz) + sqrty * sqrtz;
    sum += fac / (sqrtz * (zt * alamb));
    fac = 0.25 * fac;
    xt = 0.25 * (xt + alamb);
    yt = 0.25 * (yt + alamb);
    zt = 0.25 * (zt + alamb);
    ave = 0.2 * (xt + yt + 3 * zt);
    delx = (ave - xt) / ave;
    dely = (ave - yt) / ave;
    delz = (ave - zt) / ave;
    itC++;
    if(itC >= mIt) {
      break;
    }
  } while(Math.max(Math.max(Math.abs(delx), Math.abs(dely)), Math.abs(delz)) > minError);
  ea = delx * dely;
  eb = delz * delz;
  ec = ea - eb;
  ed = ea - 6 * eb;
  ee = ed + ec + ec;
  return 3 * sum + fac * (1 + ed * (-c1 + c5 * ed - c6 * delz * ee) + delz * (c2 * ee + delz * (-c3 * ec + delz * c4 * ea))) / (ave * ave ** 0.5);
}

function jacElliptic(u, emc) {
  let ca = 0.0003;
  let a, b, c, d;
  let em = [8];
  let en = [8];
  let bo, i, ii, I;
  let sn, cn, dn;

  if(emc != 0.0) {
    bo = 0;
    if(emc < 0) {
      bo = 1;
    }

    if(bo != 0) {
      d = 1 - emc;
      emc = -emc / d;
      d = Math.sqrt(d);
      u = d * u;
    }
    a = 1;
    dn = 1;

    for(i = 0; i < 8; i++) {
      I = i;
      em[i] = a;
      emc = Math.sqrt(emc);
      en[i] = emc;
      c = 0.5 * (a + emc);

      if(Math.abs(a - emc) <= ca * a) {
        u = c * u;
        sn = Math.sin(u);
        cn = Math.cos(u);

        if(sn == 0) {
          if(bo != 0) {
            a = dn;
            dn = cn;
            cn = a;
            sn = sn / d;
          }
          break;
        }

        a = cn / sn;
        c = a * c;
        for(ii = 1; ii >= 0; --ii) {
          b = em[ii];
          a = c * a;
          c = dn * c;
          dn = (en[ii] + a) / (b + a);
          a = c / b;
        }

        a = 1 / Math.sqrt(c * c + 1);
        if(sn < 0) {
          sn = -a;
        } else {
          sn = a;
        }
        cn = c * sn
        break;
      }
      emc = a * emc;
      a = c;
    }
  } else {
    cn = 1 / Math.cosh(u);
    dn = cn;
    sn = Math.tanh(u);
  }
  return {
    s: sn,
    c: cn,
    d: dn
  };
}

function jacobiAm(u, x, k) {
  let a = new Array(31),
    g = new Array(31),
    c = new Array(31);

  if(k == 1) {
    return 2 * Math.atan(Math.exp(u)) - Math.PI * 2;
  }

  a[0] = 1.0;
  g[0] = Math.sqrt(1.0 - k * k);
  c[0] = k;

  let two_n = 1;
  for(let n = 0; n < 30; n++) {
    if(Math.abs(a[n] - g[n]) < (a[n] * Math.EPSILON)) break;
    two_n += two_n;
    a[n + 1] = 0.5 * (a[n] + g[n]);
    g[n + 1] = Math.sqrt(a[n] * g[n]);
    c[n + 1] = 0.5 * (a[n] - g[n]);
  }
  let phi = two_n * a[n] * u;

  for(let n = 30; n > 0; n--) {
    phi = 0.5 * (phi + Math.asin(c[n] * Math.sin(phi) / a[n]));
  }
  return phi;
}

function jacobiAmA(u, x) {
  if(x == 0) {
    return u;
  }
  jacobiAm(u, x, Math.abs(x));
}

function jacobiAmM(u, x) {
  if(x == 0) {
    return u;
  }
  return jacobiAm(u, x, Math.sqrt(Math.abs(x)));
}

function jacobiAmK(u, x) {
  if(x == 0) {
    return u;
  }
  return jacobiAm(u, x, Math.sin(Math.abs(x)));
}


function addPoly(a, b) {
  if(a.length < b.length) {
    [a, b] = [b, a];
  }
  const max = a.length;
  const min = b.length;
  const result = new Array(max);
  let i = 0;
  for(; i < min; i++) {
    result[i] = a[i] + b[i];
  }
  for(; i < max; i++) {
    result[i] = a[i];
  }
  return result;
}

function multiplyPoly(a, b) {
  const al = a.length;
  if(al === 0) {
    return [];
  }
  const bl = b.length;
  if(bl === 0) {
    return [];
  }
  const size = al + bl - 2;
  const result = new Array(size).fill(0);
  let ai, bi, ri;
  let are, aim, bre, bim;
  for(ai = 0; ai < al; ai += 2) {
    are = a[ai];
    aim = a[ai + 1];
    for(bi = 0; bi < bl; bi += 2) {
      bre = b[bi];
      bim = b[bi + 1];
      ri = ai + bi;
      result[ri] += are * bre - aim * bim;
      result[ri + 1] += are * bim + bre * aim;
    }
  }
  return result;
}

function multiplyMatrices(a, b) {
  const length = a.length;
  const size = Math.sqrt(length);
  const result = new Array(length);
  for(let r = 0; r < size; r++) {
    for(let c = 0; c < size; c++) {
      let x = [];
      const offset = r * size;
      for(let i = 0; i < size; i++) {
        x = addPoly(x, multiplyPoly(a[offset + i], b[i * size + c]));
      }
      result[offset + c] = x;
    }
  }
  return result;
}

function findRoots(poly) {
  const epsilon = Number.EPSILON,
    negativeEpsilon = -Number.EPSILON;
  const length = poly.length;
  const size = length / 2;
  const roots = size - 1;
  const real = new Array(size);
  const im = new Array(size);
  for(let i = 0, j = 0; i < size; i++, j += 2) {
    real[i] = poly[j];
    im[i] = poly[j + 1];
  }
  let rc = real[roots],
    ic = im[roots],
    m = rc * rc + ic * ic;
  rc /= m;
  ic /= -m;
  let c1, c2, c3, dc = ic - rc,
    sc = rc + ic;
  for(let i = 0; i < roots; ++i) {
    c1 = rc * (real[i] + im[i]);
    c2 = real[i] * dc;
    c3 = im[i] * sc;
    real[i] = c1 - c3;
    im[i] = c1 + c2;
  }
  real[roots] = 1.0;
  im[roots] = 0.0;
  const zr = new Array(roots);
  const zi = new Array(roots).fill(0);
  for(let i = 0; i < roots; ++i) {
    zr[i] = i / 10;
  }

  let j, k, a, b, na, nb, pa, pb, qa, qb, k1, k2, k3, s1, s2, t;
  for(let i = 0; i < 1000; ++i) {
    let foundAll = true;
    for(j = 0; j < roots; ++j) {
      pa = zr[j];
      pb = zi[j];

      a = 1.0;
      b = 0.0;
      for(k = 0; k < roots; ++k) {
        if(k === j) {
          continue;
        }
        qa = pa - zr[k];
        qb = pb - zi[k];
        if(qa < epsilon && qa > negativeEpsilon && qb < epsilon && qb > negativeEpsilon) {
          continue;
        }
        k1 = qa * (a + b);
        k2 = a * (qb - qa);
        k3 = b * (qa + qb);
        a = k1 - k3;
        b = k1 + k2;
      }

      na = real[roots];
      nb = im[roots];
      s1 = pb - pa;
      s2 = pa + pb;
      for(k = size - 2; k >= 0; --k) {
        k1 = pa * (na + nb);
        k2 = na * s1;
        k3 = nb * s2;
        na = k1 - k3 + real[k];
        nb = k1 + k2 + im[k];
      }

      if(a > epsilon || a < negativeEpsilon || b > epsilon || b < negativeEpsilon) {
        k1 = a * a + b * b;
        a /= k1;
        b /= -k1;
      } else {
        a = 1.0;
        b = 0.0;
      }

      k1 = na * (a + b);
      k2 = a * (nb - na);
      k3 = b * (na + nb);

      qa = k1 - k3;
      qb = k1 + k2;

      zr[j] = pa - qa;
      zi[j] = pb - qb;

      if(qa > epsilon || qa < negativeEpsilon || qb > epsilon || qb < negativeEpsilon) {
        foundAll = false;
      }
    }

    if(foundAll) {
      break;
    }
  }

  return [zr, zi];
}

function zeta(x) {
  function sum(fn, lstRange) {
    return lstRange.reduce(
      function(lngSum, x) {
        return lngSum + fn(x);
      }, 0
    );
  }

  function range(m, n) {
    return Array.apply(null, Array(n - m + 1)).map(function(x, i) {
      return m + i;
    });
  }


  return sum(
    function(x) {
      return 1 / (x * x);
    },
    range(1, 1000)
  );
}

var factorialArray = [200];

function dynamicFactorial(n) {
  n = n < 0 ? -n : n;
  if(n == 0 || n == 1) return 1
  if(factorialArray[n] != undefined) return factorialArray[n]
  let temp = dynamicFactorial(n - 1) * n;
  factorialArray[n] = temp
  return factorialArray[n]
}

function factorial(n) {
  n = n < 0 ? -n : n
  let result = 1;
  while(n > 1) {
    result *= n;
    n--;
  }
  return result;
}

const pochhammerMap = new Map();

function dynamicPochhammer(q, n) {
  if(n == 0) return 1
  if(!pochhammerMap.has(q)) {
    pochhammerMap.set(q, [200])
  }
  if(pochhammerMap.get(q)[n] != undefined) return pochhammerMap.get(q)[n]
  let temp = dynamicPochhammer(q, n - 1) * (q + n - 1)
  pochhammerMap.get(q)[n] = temp;
  return temp
}

function pochhammer(q, n) {
  let y = 1;

  for(let k = 0; k <= n - 1; k++) {
    y *= q + k;
  }

  return y;
}

function pochhammerFalling(q, n) {
  let y = 1;

  for(let k = 0; k <= n - 1; k++) {
    y *= q - k;
  }

  return y;
}

function gammaLanczos(z) {
  let p;
  let i, y, t, x;

  p = [
    676.5203681218851,
    -1259.1392167224028,
    771.32342877765313,
    -176.61502916214059,
    12.507343278686905,
    -0.13857109526572012,
    9.9843695780195716e-6,
    1.5056327351493116e-7
  ];

  if(z < 0.5) {
    y = Math.PI / (Math.sin(Math.PI * z) * gammaLanczos(1 - z));
  } else {
    z = z - 1;
    x = 0.99999999999980993;
    for(i = 0; i < p.length; i++) {
      x = x + p[i] / (z + i + 1);
    }
    t = z + p.length - 0.5;
    y = Math.sqrt(2 * Math.PI) * t ** (z + 0.5) * Math.exp(-t) * x;
  }

  return y;
}

function tpow(z, c) {
  let a = Math.atan2(z.im, z.re);
  let lnr = 0.5 * Math.log(dot(z, z));
  let m = Math.exp(c.re * lnr - c.im * a);
  let ang = c.re * a + c.im * lnr + 2 * Math.PI;
  return multScalar(C(Math.cos(ang), Math.sin(ang)), m)
}

const gammaMap = new Map();

function dynamicGamma(z) {
  if(gammaMap.has(z)) return gammaMap.get(z)
  if((z.re < 0) && (z.re == Math.floor(z.re)) && (z.im == 0)) return C(1 / 0, 1 / 0)
  if((z.re == Math.floor(z.re)) && (z.im == 0)) return C(dynamicFactorial(z.re), 0)

  let prod = div(C(1, 0), z);
  for(let i = 1; i < 40; i++) {
    prod = mult(prod, tpow(C(1 + 1 / i, 0), z));
    prod = div(prod, addScalar(div(z, C(i, 0)), 1))
  }
  gammaMap.set(z, prod)
  return prod
}

function cGamma(z) {
  let eps = 0.1;
  let sum = C(0, 0);
  for(let i = 0; i < 10; i += eps) {
    sum = add(sum, multScalar(mult(tpow(C(i, 0), addScalar(z, -1)), exp(C(-i, 0))), eps));
  }
  return sum
}


function isInt(a) {
  return (a - Math.floor(a)) == 0;
}

function isNeg(a) {
  return a < 0;
}

function lessThan(a, b) {
  return a <= b ? a : b;
}

function greaterThan(a, b) {
  return a >= b ? a : b;
}

function binomCoeff(n, k) {
  return factorial(n) / (factorial(k) * factorial(n - k));
}

function binomCoefficient(n, k) {
  return dynamicFactorial(n) / (dynamicFactorial(k) * dynamicFactorial(n - k));
}

function round(x) {
  let f = Math.floor(x)
  let n = x > 0 ? (Math.abs(x - f) > 0.5 ? 1 : 0) : (Math.abs(x - f) > 0.5 ? -1 : 0);
  return f + n;
}

function pow(z, p) {
  if(p == 0) return C(1,0)
  return exp(multScalar(log(z), p))
}

function hypergeo2F1PowSerI(a, b, c, z) {
  let y = C(0, 0);
  let yp = C(0, 0);
  for(let k = 0; k < 40; k++) {
    yp = divScalar(multScalar(pow(z, k), pochhammer(a, k) * pochhammer(b, k)), gammaLanczos(c + k) * factorial(k));
    y = add(y, yp);
  }
  return y;
}

function powSerO(a, b, c, z) {
  let y = C(0, 0);
  let yp = C(0, 0);
  for(let k = 0; k < 40; k++) {
    yp = divScalar(multScalar(pow(z, -k), pochhammer(a, k) * pochhammer(a - c + 1, k)), factorial(k) * gammaLanczos(a - b + k + 1));
    y = add(y, yp);
  }
  return y;
}

function returnAll() {
  let result = [];
  for(let i = 0; i < arguments.length; i++) {
    result.push(arguments[i]);
  }
  return result;
}

function dmsToDec(degree, minute, second) {
  return degree + minute / 60 + second / 3600 + 1 / 1000000;
}


function cPoincare(z, p) {
  return div(add(z, p), addScalar(mult(z, conj(p)), 1))
}


function paraSize(p) {
  return (1 / (Math.sin(Math.PI * 0.5) ** 2 / Math.sin(Math.PI / p) ** 2 - 1) ** 0.5) * (Math.sin(Math.PI * 0.5) / Math.sin(Math.PI / p) - 1);
}

function square(z) {
  return mult(z,z)
}

function cWeierstrassElliptic(z, w2) {
  let w1 = C(1,0);
  let m = -100;
  let n = 100;
  let eps = 1;
  let sum = C(0,0);
  do{
    if(m == 0) {
      m = m + eps;
    }
    if(n == 0) {
      n = n - eps;
    }
    let lambda = add(multScalar(w1, m), multScalar(w2, n));
    let t = sub(div(C(1,0), square(sub(z, lambda))), div(C(1,0), square(lambda)));
    sum = add(sum, t);
    m += eps;
    n -= eps;
  } while(m < 100 && n > -100)
  return add(div(C(1,0), square(z)), sum)
}


function hypergeo2F1Coefficients(a, b, c, nomial) {
  let poly = [nomial];
  poly[0] = C(1,0);
  for(let i = 1; i < nomial; i++) {
    poly[i] = multScalar(poly[i-1], (a + i - 1) * (b + i - 1) / ((c + i - 1) * (i)));
  }
  return poly
}
var variable = 0;
function hypergeo2F1pn(z, poly) {
  let sum = C(0, 0);
  let zs = C(1,0);
  for(let i = 0; i < poly.length; i++) {
    sum = add(sum, mult(zs, poly[i]))
    zs = mult(zs, z)
  }
  variable = 10;
  return sum
}


function scHG(z, n, k, poly) {
  return mult(z, hypergeo2F1pn(pow(z, n), poly))
}


/*function derivative(z, func) {
  let h = 0.01;
  return divScalar(sub(func(addScalar(z, h)), func(z)), h)
}
function hypergeo2F1pnInv(z, poly) {
  let sum = C(0, 0);
  if(z.re == 0 &&  z.im == 0) {
    return poly[0]}
  for(let i = 0; i < poly.length; i++) {
    sum = add(sum, mult(pow(z, i), poly[i]))
    if(variable == 0) {
      console.log(sum)
    }
  }
  variable = 10;
  return sum
}*/

function nthDerivative(func, n, z) {
  if(n == 0) return func(z)
  if(n == 1) return derivative(z, func)
  else {
    function nuFunc(z) {
      return derivative(z, func)
    }
    return nthDerivative(nuFunc, n - 1, z)

  }
}

function hypergeo2F1pnInv(z, poly) {
  let sum = C(0, 0);
  if(z.re == 0 &&  z.im == 0) {
    return poly[0]}
  for(let i = 0; i < poly.length; i++) {
    sum = add(sum, mult(pow(z, i), poly[i]))
    if(variable == 0) {
      console.log(sum)
    }
  }
  variable = 10;
  return sum
}


function numericalNthDerivative(poly, n) {
  let h = 0.001;
  let sum = 0;
  for(let i = 0; i <= n; i++) {
    let pn = hypergeo2F1pn(C(i * h, 0), poly);
    pn.re = 1 / Math.pow(pn.re, n + 1) * binomCoefficient(n, i);
    pn.re *= Math.pow(-1, n + i);
    sum += pn.re;
  }
  sum = sum / Math.pow(h, n);
  return C(sum, 0)
}


function lagrangian(poly) {
  let invPoly = [poly.length];
  invPoly[0] = C(0, 0);
  for(let i = 1; i < poly.length; i++) {
    invPoly[i] = divScalar(numericalNthDerivative(poly, i - 1), dynamicFactorial(i));
  }
  return invPoly
}


function phiD(dimensions, iterations) {
  let x = 2;
  for(i = 0; i < iterations; i++) {
    x = Math.pow(1 + x, 1 / (dimensions + 1));
  }
  return x;
}

function getIrrationals(dimensions) {
  let phi = phiD(dimensions, 100);
  let ans = [];
  for(let i = 1; i <= dimensions; i++) {
    ans.push(Math.pow(1 / phi, i) % 1);
  }
  return ans;
}


/*transforms*/

function schwarzChristoffelInverseMap(n) {
  let nom = 10;
  let poly = hypergeo2F1Coefficients(1 / n, 2 / n, 1 + 1 / n, nom)
  let invPoly = lagrangian(poly);
  console.log(invPoly)
  console.log(poly)
  return z => {
    return {
      ...z,
      ...hypergeo2F1pn(z, invPoly)
    }
  }
}


function reset() {
  return z => {
    return {
      re: 0,
      im: 0,
      z: 0,
      red: 1,
      green: 1,
      blue: 1,
      alpha: 1,
    }
  }
}

function arcsinh() {
  return z => {
    return {
      ...z,
      ...multScalar(
        log(
          add(
            z,
            sqrt(addScalar(mult(z, z), 1))
          )
        ),
        2 / Math.PI)
    }
  }
}

function arctanh() {
  return z => {
    return {
      ...z,
      ...multScalar(
        log(
          div(
            addScalar(z, 1),
            addScalar(neg(z), 1)
          )
        ),
        1 / Math.PI)
    }
  }
}

function bent(real, imaginary) {
  return z => {
    let nr = z.re < 0 ? real * z.re : z.re,
      ni = z.im < 0 ? imaginary * z.im : z.im
    return {
      ...z,
      re: nr,
      im: ni
    }
  }
}

function blurCircle() {
  const [ax, ay] = getIrrationals(2);
  let atx = Math.random();
  let aty = Math.random();
  return z => {

    do {
      atx = (atx + ax) % 1;
      aty = (aty + ay) % 1;
    } while((atx - 0.5) ** 2 + (aty - 0.5) ** 2 > 0.25);

    return {
      ...z,
      re: atx * 2 - 1,
      im: aty * 2 - 1,
      z: 0
    }
  }
}

function blurGasket() {
  return z => {
    let a = Math.random() * Math.PI * 2,
      r = 1 / Math.sqrt(Math.abs(Math.random() - Math.random()));
    return {
      ...z,
      re: Math.random() - 0.5,
      im: Math.sin(a) * r
    }
  }
}

function blurGaussian(pow) {
  return z => {
    let a = gaussRnd() * Math.PI * 2,
      s = Math.sqrt(Math.abs(gaussRnd())) * pow;
    return {
      ...z,
      re: Math.cos(a) * s + z.re,
      im: Math.sin(a) * s + z.im
    }
  }
}

function blurNgon(n) {
  let points = [];
  for(let i = 0; i < n; i++) {
    points.push({
      re: Math.cos(i / n * Math.PI * 2),
      im: Math.sin(i / n * Math.PI * 2)
    });
  }

  let areaSum = 0;
  let areas = [];
  let triangles = [];
  for(let i = 2; i < points.length; i++) {
    triangles.push(blurTriangle(
      points[0], points[i - 1], points[i]));
    let area =
      Math.abs(points[0].re * points[i - 1].im +
        points[i - 1].re * points[i].im +
        points[i].re * points[0].im -
        points[0].im * points[i - 1].re -
        points[i - 1].im * points[i].re -
        points[i].im * points[0].re);
    areaSum += area;
    areas.push(area);
  }
  let len = triangles.length;

  return z => {
    let triangle = Math.random() * areaSum;
    for(let i = 0; i < len; i++) {
      if(triangle <= areas[i]) {
        return triangles[i](z);
      }
      triangle -= areas[i];
    }
  }
}

function blurSine(pow) {
  return z => {
    let a = Math.random() * 2 * Math.PI,
      u = Math.random();
    let r = (pow == 1 ? Math.acos(u * 2 - 1) : Math.acos(Math.exp(Math.log(1 - u) * pow) * 2 - 1)) / Math.PI;
    return {
      ...z,
      re: Math.cos(a) * r,
      im: Math.sin(a) * r
    }
  }
}

function blurSquare() {
  const [ax, ay] = getIrrationals(2);
  let atx = Math.random();
  let aty = Math.random();
  return z => {
    atx = (atx + ax) % 1;
    aty = (aty + ay) % 1;

    return {
      ...z,
      re: atx - 0.5,
      im: aty - 0.5,
      z: 0,
    }
  }
}

function squircleDistance(z, r) {
    if(z.im > -r && z.im < r) {
      if(z.re < -r) {
        return -r - z.re;
      }
      if(z.re > r) {
        return z.re - r;
      }
      return 0;
    }
    if(z.re > -r && z.re < r) {
      if(z.im < -r) {
        return -r - z.im;
      }
      if(z.im > r) {
        return z.im - r;
      }
      return 0;
    }

    return Math.hypot(Math.abs(z.re) - r, Math.abs(z.im) - r);
}

function blurTriangle(a, b, c) {
  let nb = sub(b, a);
  let nc = sub(c, a);
  return z => {
    let t1 = Math.random();
    let t2 = Math.random();
    if(t1 + t2 > 1) {
      t1 = 1 - t1;
      t2 = 1 - t2;
    }
    let mb = multScalar(nb, t1)
    let mc = multScalar(nc, t2);
    let f = add(add(mb, mc), a)
    return {
      ...z,
      ...f
    }
  }
}

function triangle3D(x1, y1, z1, x2, y2, z2, x3, y3, z3) {
  let t1 = Math.random();
  let t2 = Math.random();

  if(t1 + t2 > 1) {
    t1 = 1 - t1;
    t2 = 1 - t2;
  }

  return {
    re: x1 + (x2 - x1) * t1 + (x3 - x1) * t2,
    im: y1 + (y2 - y1) * t1 + (y3 - y1) * t2,
    z: z1 + (z2 - z1) * t1 + (z3 - z1) * t2,
  };
}

function bTransform(rotate, power, move, split) {
  return z => {
    let tau = 0.5 * (Math.log((z.re + 1) * (z.re + 1) + z.im * z.im) - Math.log((z.re - 1) * (z.re - 1) + z.im * z.im)) / power + move,
      sigma = Math.PI - Math.atan2(2 * z.im, 1 - dot(z, z)) + rotate;
    let sigma2 = (sigma + 2 * Math.PI * Math.floor(Math.random() * power)) / power;
    tau2 = z.re >= 0 ? tau + split : tau - split;
    let f = Math.cosh(tau2) - Math.cos(sigma2)
    return {
      ...z,
      re: Math.sinh(tau2) / f,
      im: Math.sin(sigma2) / f
    }
  }
}

function bubble() {
  return z => {
    let r = 4 / (dot(z, z) + 4)
    return {
      ...z,
      re: z.re * r,
      im: -z.im * r
    }
  }
}

function circleInv() {
  return z => {
    let r = 1 / dot(z, z)
    return {
      ...z,
      re: z.re * r,
      im: -z.im * r
    }
  }
}

function cpow(rotation, depth, divisor, spread) {
  let ang = 2 * Math.PI / divisor,
    a1 = depth < 0 ? -Math.log(-depth) * rotation : Math.log(depth) * rotation;
  let a = Math.atan2(a1, 2 * Math.PI);
  let c = Math.cos(a) * rotation * Math.cos(a) / divisor,
    d = Math.cos(a) * rotation * Math.sin(a) / divisor;
  let coeff = d == 0 ? 0 : -0.095 * spread / d;
  return z => {
    let a0 = Math.atan2(z.im, z.re);
    let a1 = a0 < 0 ? a0 + 2 * Math.PI : a0;
    let a2 = Math.cos(a1 * 0.5) < Math.random() * 2 - 1 ? a1 - 2 * Math.PI : a1;
    let a = a2 + ((Math.random() > 0.5 ? 2 : -2) * Math.PI) * Math.floor(Math.log(Math.random()) * coeff + 0.5),
      lnr2 = Math.log(dot(z, z));
    let r = Math.exp(c * 0.5 * lnr2 - d * a),
      angle = c * a + d * 0.5 * lnr2 + ang * Math.floor(Math.random() * 32767);
    return {
      ...z,
      re: Math.sin(angle) * r,
      im: Math.cos(angle) * r
    }
  }
}

function cylinder() {
  return z => {
    return {
      ...z,
      re: z.re / Math.sqrt(z.re * z.re + 1),
      im: z.im
    }
  }
}

function dc_poincareDisc(p, q, iterations, checks0, checks1, checks2) {
  let mba = C(1, 0);
  let mbb = C(0, 0);
  let a = Math.PI / p;
  let b = Math.PI / q;
  if(a + b >= 0.5 * Math.PI) {
    b = 0.5 * Math.PI - a - 0.01;
  }
  let ca = Math.cos(a),
    sa = Math.sin(a);
  let sb = Math.sin(b);
  let c1 = Math.cos(a + b);
  let s = c1 / (1 - sa * sa - sb * sb) ** 0.5;
  let sx = (s * s + 1) / (2 * s * ca);
  let sr = (sx * sx - 1) ** 0.5;
  let nfp = C(sa, -ca);
  let angle = 29.9664;
  return z => {
    let col = {
      red: z.red,
      green: z.green,
      blue: z.blue
    };
    if(dot(z, z) > 1) {
      col = {
        red: 0,
        green: 0,
        blue: 0
      };
    }
    let nz = discShift(z, mba, mbb);
    let az = addScalar(nz, 1);
    let fc = {
      re: 0,
      im: 0,
      w: 0
    };
    for(let i = 0;
      (i < iterations) && (az.re == nz.re) && (az.im == nz.im); i++) {
      az = nz;
      nz.im = Math.abs(nz.im);
      fc.re += (az.im != nz.im) ? 1 : 0;
      let t = 2 * Math.min(0, dot(nz, nfp));
      fc.im += (t < 0) ? 1 : 0;
      nz = add(nz, neg(multScalar(nfp, t)));
      nz.re = sx;
      let r2 = dot(nz, nz);
      let k = Math.max(sr * sr / r2, 1);
      fc.w += (k != 1) ? 1 : 0;
      nz = multScalar(nz, k);
      nz.re += sx;
    }
    let r = Math.pow(Math.abs(Math.hypot(z.re - sx, z.im) - sr), 0.25);
    let k = ((checks0 == 1 ? 1 : 0) * fc.re + (checks1 == 1 ? 1 : 0) * fc.im + (checks2 == 1 ? 1 : 0) * fc.w) % 2;
    let k2 = 0.5 * k + 0.5;
    col = {
      red: r * k2,
      green: r * k2,
      blue: r * k2
    };
    let fz = nz;
    fz.im *= k * 2 - 1;
    fz = mult(fz, mba);
    fz = mult(fz, C(Math.cos(angle * DEGREE), Math.sin(angle * DEGREE)));
    fz = add(fz, neg(mbb));
    fz = addScalar(multScalar(fz, 0.5), 0.5);
    return {
      ...z,
      ...fz,
      red: col.red,
      green: col.green,
      blue: col.blue
    }
  }
}

function disc() {
  return z => {
    let r = Math.sqrt(dot(z, z)) * Math.PI,
      a = Math.atan2(z.im, z.re) / Math.PI;
    return {
      ...z,
      re: Math.sin(r) * a,
      im: Math.cos(r) * a
    }
  }
}

function dragon(a, divisorB, divisorC, bc, multiplier, horizontal, vertical, radial) {
  bc *= 2;
  const offsetB = 2 * Math.cos(Math.PI / divisorB),
    offsetC = 2 * Math.cos(Math.PI / divisorC);
  let real = 0,
    im = 0;
  let preOffsets = [];
  let postOffsets = [];

  if(bc === 0) {
    real = 2 * Math.cos(Math.PI / multiplier / a);
    preOffsets[0] = 0;
    postOffsets[0] = 0;
  } else {
    let matrix = [
        [1, 0],
        [],
        [],
        [1, 0]
      ],
      n = a,
      b = true,
      i = 0,
      currentOffset = 0,
      offset;
    for(let x = 0; x < bc; x++, n += a) {
      while(n > 0) {
        n -= bc;
        matrix = multiplyMatrices([
          [],
          [1, 0],
          [-1, 0],
          [0, 0, 1, 0]
        ], matrix);
        preOffsets[i] = currentOffset;
        postOffsets[i++] = currentOffset;
      }
      offset = b ? -offsetB : -offsetC;
      matrix = multiplyMatrices([
        [],
        [1, 0],
        [-1, 0],
        [0, offset, 1, 0]
      ], matrix);
      preOffsets[i] = offset + currentOffset;
      if(currentOffset !== 0) {
        currentOffset = 0;
      } else {
        currentOffset = offset;
      }
      postOffsets[i++] = currentOffset;
      b = !b;
    }
    const trace = addPoly(matrix[0], matrix[3]);
    trace[0] -= 2 * Math.cos(Math.PI / multiplier);
    const [realParts, imParts] = findRoots(trace);
    for(let i = 0, len = realParts.length; i < len; i++) {
      if(realParts[i] > real) {
        real = realParts[i];
        im = imParts[i];
      }
    }
  }
  for(let i = 0; i < preOffsets.length; i++) {
    preOffsets[i] += im;
  };

  let current = 0,
    rotate = true,
    direct = true;
  const rotationOffset = im - offsetB;
  const length = preOffsets.length;

  const offset = offsetB + offsetC;
  return z => {
    let x = z.re,
      y = z.im;
    if(Math.random() > horizontal) {
      direct = Math.random() > 0.5;
      current = Math.floor(Math.random() * length);
      if(Math.random() > 0.5) {
        rotate = !rotate;
      }
    }
    if(Math.random() > radial) {
      rotate = !rotate;
      current = direct ? ++current % length : ((current === 0) ? length - 1 : --current);
    }
    if(direct) {
      if(rotate) {
        y += preOffsets[current] - rotationOffset;
      } else {
        x = real - x;
        y = preOffsets[current] - y;
      }
      const c = 1 / (x * x + y * y);
      x = x * c;
      y = postOffsets[current] - y * c;
    } else {
      if(rotate) {
        x = real - x;
        y = rotationOffset - y - postOffsets[current];
      } else {
        y -= postOffsets[current];
      }
      const c = 1 / (x * x + y * y);
      x = real - x * c;
      y = preOffsets[current] + y * c;
    }
    y += offset * (Math.round(Math.cos(Math.random() * Math.PI) * (1 / Math.sqrt(Math.random()) - Math.random()) * vertical / offset));
    if(Math.random() > 0.5) {
      return {
        ...z,
        re: real - x,
        im: rotationOffset - y
      }
      rotate = false;
    } else {
      return {
        ...z,
        re: x,
        im: y
      }
      rotate = true;
    }
  }
}

function ePush(push, _rotation) {
  rotation = _rotation * DEGREE;
  return z => {
    let tmp = dot(z, z) + 1;
    let tmp2 = 2 * z.re;
    let xm = ((tmp + tmp2) ** 0.5 + (tmp - tmp2) ** 0.5) * 0.5;
    let xmax = (xm < 1) ? 1 : xm;
    let nt = z.re / xmax;
    let t = (nt > 1) ? 1 : (nt < -1) ? -1 : nt;

    let nu = ((z.im < 0) ? -Math.acos(t) : Math.acos(t)) + rotation;
    let mu = Math.acosh(xmax) + push;
    return {
      ...z,
      re: Math.cosh(mu) * Math.cos(nu),
      im: Math.sinh(mu) * Math.sin(nu)
    }
  }
}

function eRotate(rotationAngle) {
  rotationAngle *= DEGREE;
  return z => {
    let tmp = dot(z, z) + 1;
    let tmp2 = 2 * z.re;
    let xm = ((tmp + tmp2) ** 0.5 + (tmp - tmp2) ** 0.5) * 0.5;
    xm = xm < 1 ? 1 : xm;
    let t = z.re / xm;
    t = t > 1 ? 1 : t < -1 ? -1 : t;
    let nu = z.im < 0 ? -Math.acos(t) : Math.acos(t);
    nu = (nu + rotationAngle + Math.PI) % (2 * Math.PI) - Math.PI;
    return {
      ...z,
      re: xm * Math.cos(nu),
      im: (xm * xm - 1) ** 0.5 * Math.sin(nu)
    }
  }
}


function flipX() {
  return z => {
    return {
      ...z,
      re: z.im > 0 ? -z.re : z.re
    }
  }
}

function flipY() {
  return z => {
    return {
      ...z,
      im: z.re > 0 ? -z.im : z.im
    }
  }
}

function hypershape(n) {
  let alpha = 2 * Math.PI / n;
  let alphacoeff = Math.tan(alpha * 0.5) * 2,
    beta = Math.SQRT2 * Math.cos(alpha * 0.5);
  return z => {
    let da = (Math.atan2(z.im, z.re) + Math.PI) / alpha,
      rad = Math.sqrt(dot(z, z));
    let za0 = Math.floor(da);
    let xa0 = da - za0,
      xa1, za, si;
    if(xa0 > 0.5) {
      xa1 = 1 - xa0;
      za = za0 + 1;
      si = -1;
    } else {
      xa1 = xa0;
      za = za0;
      si = 1;
    }
    let xa = Math.atan(xa1 * alphacoeff) / alpha;
    let co = beta / Math.cos(xa * alpha),
      ang = (za + si * xa) * alpha - Math.PI;

    let p = multScalar({
      re: Math.cos(ang),
      im: Math.sin(ang)
    }, co * rad);
    let r = dot(p, p);

    return {
      ...z,
      ...(divScalar(multScalar(p, Math.SQRT2 - Math.sqrt(2 - r)), r))
    }
  }
}

function hypershift(p) {
  return z => {
    return {
      ...z,
      ...div(add(z, p), addScalar(mult(conj(p), z), 1))
    }
  }
}

function hypertile3(p, q, r, shift) {
  let o = Math.acosh((Math.cos(Math.PI / p) + Math.cos(Math.PI / q) * Math.cos(Math.PI / r)) / (Math.sin(Math.PI / q) * Math.sin(Math.PI / r)));
  let a = Math.asinh(Math.sin(Math.PI / q) / Math.sin(Math.PI / p) * Math.sinh(o));
  let b = Math.asinh(Math.sin(Math.PI / r) * Math.sinh(a)),
    c = Math.asinh(Math.sin(Math.PI / r) * Math.sinh(a) / Math.sin(Math.PI / q));
  let h = Math.tanh(b / 2),
    b1 = Math.tanh(Math.acosh(Math.cosh(c) / Math.cosh(b)) / 2),
    b2 = Math.tanh(Math.acosh(Math.cosh(a) / Math.cosh(b)) / 2);

  let c1 = C(h, 0),
    c2 = C(0, b1),
    c3 = C(0, b2);

  return z => {
    let rot1 = 360 / p * Math.floor(Math.random() * p) * DEGREE,
      rot2 = 360 / q * Math.floor(Math.random() * q) * DEGREE,
      rot3 = 360 / r * Math.floor(Math.random() * r) * DEGREE;
    let r01 = C(Math.cos(rot1), Math.sin(rot1)),
      r02 = C(Math.cos(rot2), Math.sin(rot2)),
      r03 = C(Math.cos(rot3), Math.sin(rot3));
    let n, pfr, z0;
    if(shift < 0.25) {
      n = 0;
      pfr = 0;
      z0 = z;
    } else if(shift < 0.5) {
      n = p;
      pfr = rot1;
      z0 = div(add(z, c1), addScalar(mult(c1, z), 1));
    } else if(shift < 0.75) {
      n = q;
      pfr = rot2;
      z0 = div(add(z, neg(c2)), addScalar(mult(c2, z), 1));
    } else {
      n = r;
      pfr = rot3;
      z0 = div(add(z, c3), addScalar(mult(neg(c3), z), 1));
    }

    let m0 = div(add(z0, neg(c1)), addScalar(mult(neg(c1), z0), 1));
    m0 = mult(r01, m0);
    m0 = div(add(m0, c1), addScalar(mult(c1, m0), 1));
    let m1 = div(add(z0, c2), addScalar(mult(neg(c2), z0), 1));
    m1 = mult(r02, m1);
    m1 = div(add(m1, neg(c2)), addScalar(mult(c2, m1), 1));
    let m2 = div(add(z0, neg(c3)), addScalar(mult(c3, z0), 1));
    m2 = mult(r03, m2);
    m2 = div(add(m2, c3), addScalar(mult(neg(c3), m2), 1));

    let fr = Math.floor(Math.random() * n) * pfr,
      rnd = Math.random(),
      f3, f0, f;

    if(rnd < 1 / 3) {
      f3 = m0
    } else if(rnd < 2 / 3) {
      f3 = m1
    } else {
      f3 = m2
    }

    if(shift < 0.25) {
      f0 = f3
    } else if(shift < 0.5) {
      f0 = div(add(f3, neg(c1)), addScalar(mult(neg(c1), f3), 1))
    } else if(shift < 0.75) {
      f0 = div(add(f3, c2), addScalar(mult(neg(c2), f3), 1))
    } else {
      f0 = div(add(f3, neg(c3)), addScalar(mult(c3, f3), 1))
    }

    if(shift < 0.25) {
      f = f0
    } else {
      f = mult(C(Math.cos(fr), Math.sin(fr)), f0)
    }

    return {
      ...z,
      ...f
    }
  }
}

function jac_cd(k) {
  return z => {
    let jx = jacElliptic(z.re, k);
    let jy = jacElliptic(z.im, 1 - k);

    let numxc = jx.c * jy.c;
    let numyc = jx.d * jx.s * jy.d * jy.s;

    let numxd = jx.d * jy.c * jy.d;
    let numyd = jx.c * jx.s * jy.s * k;

    let denom = 1 / (jx.s * jx.s * jy.s * jy.s * k + jy.c * jy.c);

    let cn = {re: numxc * denom, im: numyc * denom}
    let dn = {re: numxd * denom, im: numyd * denom}

    return {
      ...z,
      ...div(cn, dn)
    }
  }
}

function jac_cn(k) {
  return z => {
    let jx = jacElliptic(z.re, k);
    let jy = jacElliptic(z.im, 1 - k);

    let numx = jx.c * jy.c;
    let numy = jx.d * jx.s * jy.d * jy.s;

    let denom = 1 / (jx.s * jx.s * jy.s * jy.s * k + jy.c * jy.c);

    return {
      ...z,
      re: denom * numx,
      im: denom * numy
    }
  }
}

function jac_dn(k) {
  return z => {
    let jx = jacElliptic(z.re, k);
    let jy = jacElliptic(z.im, 1 - k);

    let numx = jx.d * jy.c * jy.d;
    let numy = jx.c * jx.s * jy.s * k;

    let denom = 1 / (jx.s * jx.s * jy.s * jy.s * k + jy.c * jy.c);

    return {
      ...z,
      re: denom * numx,
      im: denom * numy
    }
  }
}

function jac_elk(k) {
  let ephi, epsi, sinA, cosA, phi, psi;

  return z => {
    phi = z.re;
    psi = z.im;

    let cotphi2 = ssqr(scot(phi));

    let b = -(cotphi2 + k * ssqr(Math.sinh(psi) / (Number.EPSILON + Math.sin(phi))) - 1 + k);
    b = b * 0.5;

    let c = -(1 - k) * cotphi2;
    c = ssqrt2(ssqr(b) - c);

    let x1 = Math.max(-b + c, -b - c);
    let mu = ssqrt2((x1 / nonz(cotphi2) - 1) / nonz(k));
    mu = sign(psi) * mu;
    let lambda = sign(phi) * ssqrt2(x1);

    sinA = sign(lambda) / Math.sqrt(ssqr(lambda) + 1);
    cosA = lambda * sinA;
    ephi = sinA * carlsonRF(ssqr(cosA), 1 - k * ssqr(sinA), 1);

    cosA = 1 / Math.sqrt(ssqr(mu) + 1);
    sinA = mu * cosA;
    epsi = sinA * carlsonRF(ssqr(cosA), 1 - (1 - k) * ssqr(sinA), 1);
    return {
      ...z,
      re: ephi,
      im: epsi
    }
  }
}


function jac_sn(k) {
  return z => {
    let jx = jacElliptic(z.re, k);
    let jy = jacElliptic(z.im, 1 - k);

    let numx = jx.s * jy.d;
    let numy = jx.c * jx.d * jy.c * jy.s;

    let denom = 1 / (jx.s * jx.s * jy.s * jy.s * k + jy.c * jy.c);

    return {
      ...z,
      re: denom * numx,
      im: denom * numy
    }
  }
}

function julian(pow, dist) {
  return z => {
    let a = (Math.atan2(z.im, z.re) + Math.floor(pow * Math.random()) * Math.PI * 2.0) / pow,
      r = Math.pow(dot(z, z), dist / pow * 0.5);
    return {
      ...z,
      re: Math.cos(a) * r,
      im: Math.sin(a) * r
    }
  }
}

function juliaq(pow, div) {
  let ip = div / pow,
    ip2p = (2 * Math.PI) / pow;
  return z => {
    let ang = Math.atan2(z.im, z.re) * ip + Math.floor(32767 * Math.random()) * ip2p;
    let cosa = Math.cos(ang),
      sina = Math.sin(ang),
      r = Math.pow(dot(z, z), 0.5 * ip);
    return {
      ...z,
      re: cosa * r,
      im: sina * r
    }
  }
}

function juliascope(pow, dist) {
  let ndist = dist / pow * 0.5;
  return z => {
    let root = Math.floor(pow * Math.random());
    let rootd2 = Math.floor(root * 0.5);
    let roots = root - rootd2 * 2 == 1 ? -1 : 1;
    let a = (Math.atan2(z.im, z.re) * roots + root * 2 * Math.PI) / pow,
      r = Math.pow(dot(z, z), ndist);
    return {
      ...z,
      re: Math.cos(a) * r,
      im: Math.sin(a) * r
    }
  }
}

function mobius(a, b, c, d) {
  return z => {
    return {
      ...z,
      ...div(add(mult(a, z), b), add(mult(c, z), d))
    }
  }
}

function multiMobius(a, b, c, d, iterations) {
  let mo = mobius(a, b, c, d);
  return z => {
    let nz = z;
    for(let i = 0; i < iterations; i++) {
      nz = mo(nz);
    }
    return {
      ...z,
      ...nz
    }
  }
}

function murl2(c, pow) {
  return z => {
    let angle = Math.atan2(z.im, z.re) * pow;
    let cosa = Math.cos(angle),
      sina = Math.sin(angle),
      r = c * Math.pow(z.re * z.re + z.im * z.im, 0.5 * pow);
    let real = r * cosa + 1,
      imag = r * sina;
    let r1 = Math.pow(real * real + imag * imag, 1.0 / pow),
      angle1 = Math.atan2(imag, real) * 2 / pow;
    let cosa1 = Math.cos(angle1),
      sina1 = Math.sin(angle1);
    let re2 = r1 * cosa1,
      im2 = r1 * sina1,
      vp = c != -1.0 ? Math.pow(c + 1, 2.0 / pow) : 0.0;
    let r2 = vp / (r1 * r1);
    return {
      ...z,
      re: (z.re * re2 + z.im * im2) * r2,
      im: (z.im * re2 - z.re * im2) * r2
    };
  }
}

function nSplit(n, split, wedge) {
  let regSize = Math.PI * 2 / n;

  return z => {
    let rad = Math.hypot(z.re, z.im);

    let theta = Math.atan2(z.im, z.re) + Math.PI;
    let region = Math.floor(theta / (Math.PI * 2) * n);
    let mTheta = theta - regSize * region;

    let newTheta = region * regSize + regSize / 2 + (mTheta - regSize / 2) * wedge;

    let cosA = Math.cos((region + 0.5) / n * 2 * Math.PI) * split;
    let sinA = Math.sin((region + 0.5) / n * 2 * Math.PI) * split;

    return {
      ...z,
      re: Math.cos(newTheta) * rad + cosA,
      im: Math.sin(newTheta) * rad + sinA
    }
  }
}

function pdj(a, b, c, d) {
  return z => {
    let ny1 = Math.sin(a * z.im);
    let nx1 = Math.cos(b * z.re);
    let nx2 = Math.sin(c * z.re);
    let ny2 = Math.cos(d * z.im);
    return {
      ...z,
      re: ny1 - nx1,
      im: nx2 - ny2
    }
  }
}

function pointSymmetry(centerX, centerY, order) {
  return z => {
    let idr = Math.floor(Math.random() * order),
      dx = z.re - centerX,
      dy = z.im - centerY,
      da = (2 * Math.PI) / order;
    let angle = idr * da;
    let cosa = Math.cos(angle),
      sina = Math.sin(angle)
    return {
      ...z,
      re: centerX + dx * cosa + dy * sina,
      im: centerY + dy * cosa - dx * sina
    }
  }
}

function rotate(theta) {
  const rad = DEGREE * theta;
  const sinTheta = -Math.sin(rad);
  const cosTheta = -Math.cos(rad);
  return z => {
    return {
      ...z,
      re: -(cosTheta * z.re + sinTheta * z.im),
      im: sinTheta * z.re - cosTheta * z.im
    }
  }
}

function identity() {
  return z => {
    return z;
  }
}

function scale(factor) {
  return function(z) {
    return {
      ...z,
      re: z.re * factor,
      im: z.im * factor
    }
  }
}

function scale2(real, imaginary) {
  return function(z) {
    return {
      ...z,
      re: z.re * real,
      im: z.im * imaginary
    }
  }
}

function schwarzChristoffelMap(n, k) {
  let nom = 500;
  let poly = hypergeo2F1Coefficients(1 / n, 2 / n, 1 + 1 / n, nom)
  console.log(poly)
  return z => {
    return {
      ...z,
      ...scHG(z, n, k, poly)
    }
  }
}

function schwarzTriangle(_alph, _bet, _gam) {
  let alph = _alph * DEGREE;
  let bet = _bet * DEGREE;
  let gam = _gam * DEGREE;
  let a = (1 - alph - bet - gam) / 2;
  let b = (1 - alph + bet - gam) / 2;
  let c = 1 - alph;
  let na = a - c + 1;
  let nb = b - c + 1;
  let nc = 2 - c;

  return z => {
    let nz = pow(z, alph);
    let hg = hypergeo2F1(a, b, c, z);
    let nhg = hypergeo2F1(na, nb, nc, z);

    return {
      ...z,
      ...mult(nz, div(nhg, hg))
    }
  }
}

function sinusoidal() {
  return z => {
    return {
      ...z,
      re: Math.sin(z.re),
      im: Math.sin(z.im)
    }
  }
}

function skew(real, imaginary) {
  return z => {
    return {
      ...z,
      re: z.re + real * z.im,
      im: imaginary * z.re + z.im
    }
  }
}

function smartcrop(power, radius, roundstr, roundwidth, distortion, cropmode) {
  let mode = (power > 0) == (radius > 0) ? 1 : 0,
    pow = Math.abs(power);
  let alpha = Math.PI * 2 / pow;
  let roundcoeff = roundstr / Math.sin(alpha * 0.5) / pow * 2,
    wradius = Math.abs(radius),
    radial, wpower, walpha, wcoeff;
  if(pow < 2) {
    radial = 1;
    wpower = pow * Math.PI;
    walpha = 0;
    wcoeff = 0;
  } else {
    radial = 0;
    wpower = pow;
    walpha = alpha;
    wcoeff = roundcoeff;
  }
  return z => {
    let ang = Math.atan2(z.im, z.re),
      rad = Math.sqrt(dot(z, z));
    let wedge, xang0, xang1, xang, coeff0, coeff1, coeff, xr, angle, wwidth, rdc;
    if(radial == 1) {
      xang0 = ang / (2 * Math.PI) + 1;
      xang = (xang0 - Math.floor(xang0)) * 2 * Math.PI;
      angle = Math.floor(Math.random() * 2) != 0 ? wpower : 0;
      if((xang > wpower) == (mode != 1)) {
        return {
          ...z,
          re: Math.cos(angle) * rad,
          im: Math.sin(angle) * rad
        }
      } else {
        return {
          ...z
        }
      }
    } else {
      xang0 = (ang + Math.PI) / walpha;
      xang1 = xang0 - Math.floor(xang0);
      xang = xang1 < 0.5 ? xang1 : 1 - xang1;
      coeff0 = 1 / Math.cos(xang * walpha);
      wwidth = roundwidth != 1 ? Math.exp(Math.log(xang * 2) * roundwidth) * wcoeff : xang * 2 * wcoeff;
      coeff1 = distortion == 0 ? 1 : roundstr != 0 ? Math.abs((1 - wwidth) * coeff0 + wwidth) : coeff0;
      coeff = distortion != 1 ? Math.exp(Math.log(coeff1) * distortion) : coeff1;
      xr = coeff * wradius;
      rdc = cropmode == -1 ? rad : xr;
      f = (rad > xr) == (mode == 1) ? cropmode != 0 ? multScalar(C(Math.cos(ang), Math.sin(ang)), rdc) : C(0, 0) : z;
      return {
        ...z,
        ...f
      }
    }
  }
}

function smartshape(power, roundstr, roundwidth, distortion, compensation) {
  let pow = Math.max(power, 2);
  let alpha = Math.PI * 2 / pow;
  let alphacoeff = Math.tan(alpha * 0.5) * 2,
    roundcoeff = roundstr / Math.sin(alpha * 0.5) / pow * 2,
    comp = compensation <= 0 ? 0 : 1;

  return z => {
    let dang = (Math.atan2(z.im, z.re) + Math.PI) / alpha,
      rad = Math.sqrt(dot(z, z));
    let zang1 = Math.floor(dang);
    let xang1 = dang - zang1,
      xang2, zang, sign, xang;

    if(xang1 > 0.5) {
      xang2 = 1 - xang1;
      zang = zang1 + 1;
      sign = -1;
    } else {
      xang2 = xang1;
      zang = zang1;
      sign = 1;
    }
    if(comp == 1 && distortion >= 1) {
      xang = Math.atan(xang2 * alphacoeff) / alpha;
    } else {
      xang = xang2;
    }
    let coeff0 = 1 / Math.cos(xang * alpha),
      wwidth = roundwidth != 1 ? Math.exp(Math.log(xang * 2) * roundwidth) * roundcoeff : xang * 2 * roundcoeff;
    let coeff1 = distortion == 0 ? 1 : roundstr != 0 ? Math.abs((1 - wwidth) * coeff0 + wwidth) : coeff0;
    let coeff = distortion != 1 ? Math.exp(Math.log(coeff1) * distortion) : coeff1;
    let ang = (zang + sign * xang) * alpha - Math.PI;

    return {
      ...z,
      re: Math.cos(ang) * coeff * rad,
      im: Math.sin(ang) * coeff * rad
    }
  }
}

function splits(real, imaginary) {
  return z => {
    let xoff = z.re > 0 ? real : -real;
    let yoff = z.im > 0 ? imaginary : -imaginary;
    return {
      ...z,
      re: z.re + xoff,
      im: z.im + yoff
    }
  }
}

function tileHelp() {
  return z => {
    let x = z.re;
    let val = Math.cos((x > 0 ? x - Math.floor(x) : x + Math.floor(-x)) * Math.PI),
      fpx;
    if(val < Math.random() * 2 - 1) {
      fpx = x > 0 ? -1 : 1
    } else {
      fpx = 0
    }
    return {
      ...z,
      re: z.re + fpx,
      im: z.im
    }
  }
}

function tileLog(spread) {
  return z => {
    spread = spread == 0 ? 1 / Math.sqrt(Math.abs(Math.random() - Math.random())) : spread;
    return {
      ...z,
      re: z.re + Math.floor(Math.log(Math.random()) * (Math.random() < 0.5 ? spread : -spread) + 0.5),
      im: z.im
    }
  }
}

function translate(real, imaginary) {
  return z => {
    return {
      ...z,
      re: z.re + real,
      im: z.im + imaginary
    }
  }
}

function trigCosh() {
  return z => {
    return {
      ...z,
      ...cosh(z)
    }
  }
}

function trigExp() {
  return z => {
    return {
      ...z,
      ...exp(z)
    }
  }
}

function trigLog() {
  return z => {
    return {
      ...z,
      ...log(z)
    }
  }
}

function trigSinh() {
  return z => {
    return {
      ...z,
      ...sinh(z)
    }
  }
}

function trigTanh() {
  return z => {
    return {
      ...z,
      ...tanh(z)
    }
  }
}

function unbubble() {
  return z => {
    let r = dot(z, z);
    let b = (Math.SQRT2 - Math.sqrt(2 - r)) / r;
    return {
      ...z,
      re: z.re * b,
      im: z.im * b
    }
  }
}

function weierstrassElliptic(w2) {
  return z => {
    return {
      ...z,
      ...cWeierstrassElliptic(z, w2)
    }
  }
}

/* colors */

function color(color) {
  let col = color;
  if(color.hasOwnProperty('h')) {
    col = hslToRgb(color.h, color.s, color.l);
  } else if(!color.hasOwnProperty('red')) {
    col = colorRGB(...arguments);
  }
  return z => {
    return {
      ...z,
      red: col.red,
      green: col.green,
      blue: col.blue
    }
  }
}

function lerp(colorA, colorB, weight = 0.5) {
  weight = Math.max(0, Math.min(1, weight));
  if(colorA.hasOwnProperty('h')) {
    if(colorB.hasOwnProperty('h')) {
      return hslToRgb(
        Math.abs((colorA.h % 1) - (colorB.h % 1)) > 0.5 ?
        (((colorA.h + (colorA.h % 1 > colorB.h ? 0 : 1)) * (1 - weight) + (colorB.h + (colorA.h % 1 > colorB.h ? 1 : 0)) * weight + 1) + 1) % 1 :
        colorA.h * (1 - weight) + colorB.h * weight,
        colorA.s * (1 - weight) + colorB.s * weight,
        colorA.l * (1 - weight) + colorB.l * weight);
    } else {
      colorB = rgbToHsl(colorB.red, colorB.green, colorB.blue);
      return hslToRgb(
        Math.abs((colorA.h % 1) - (colorB.h % 1)) > 0.5 ?
        (((colorA.h + (colorA.h % 1 > colorB.h ? 0 : 1)) * (1 - weight) + (colorB.h + (colorA.h % 1 > colorB.h ? 1 : 0)) * weight + 1) + 1) % 1 :
        colorA.h * (1 - weight) + colorB.h * weight,
        colorA.s * (1 - weight) + colorB.s * weight,
        colorA.l * (1 - weight) + colorB.l * weight);
    }
  } else {
    if(colorB.hasOwnProperty('h')) {
      colorA = rgbToHsl(colorA.red, colorA.green, colorA.blue);
      return hslToRgb(
        Math.abs((colorA.h % 1) - (colorB.h % 1)) > 0.5 ?
        (((colorA.h + (colorA.h % 1 > colorB.h ? 0 : 1)) * (1 - weight) + (colorB.h + (colorA.h % 1 > colorB.h ? 1 : 0)) * weight + 1) + 1) % 1 :
        colorA.h * (1 - weight) + colorB.h * weight,
        colorA.s * (1 - weight) + colorB.s * weight,
        colorA.l * (1 - weight) + colorB.l * weight);
    } else {
      return lerpRGB(colorB.red, colorB.green, colorB.blue, weight)(colorA);
    }
  }
}

function lerpColor(color, weight) {
  return z => {
    return {
      ...z,
      ...lerp(z, color, weight)
    }
  }
}

function colorRGB(r, g, b) {
  return {
    red: r,
    green: g,
    blue: b
  }
}

function colorHSL(h, s, l) {
  return {
    h: h,
    s: s,
    l: l
  };
}

function hue2rgb(p, q, t) {
  if(t < 0) t += 1;
  if(t > 1) t -= 1;
  if(t < 1 / 6) return p + (q - p) * 6 * t;
  if(t < 1 / 2) return q;
  if(t < 2 / 3) return p + (q - p) * (2 / 3 - t) * 6;
  return p;
}

function hslToRgb(h, s, l) {
  let r, g, b;

  if(s == 0) {
    r = g = b = l; // achromatic
  } else {
    let q = l < 0.5 ? l * (1 + s) : l + s - l * s;
    let p = 2 * l - q;

    r = hue2rgb(p, q, h + 1 / 3);
    g = hue2rgb(p, q, h);
    b = hue2rgb(p, q, h - 1 / 3);
  }

  return {
    red: r,
    green: g,
    blue: b
  };
}

function rgbToHsl(r, g, b) {
  let max = Math.max(r, g, b),
    min = Math.min(r, g, b);
  let h, s, l = (max + min) / 2;

  if(max == min) {
    h = s = 0; // achromatic
  } else {
    let d = max - min;
    s = l > 0.5 ? d / (2 - max - min) : d / (max + min);

    switch (max) {
      case r:
        h = (g - b) / d + (g < b ? 6 : 0);
        break;
      case g:
        h = (b - r) / d + 2;
        break;
      case b:
        h = (r - g) / d + 4;
        break;
    }

    h /= 6;
  }

  return {
    h: h,
    s: s,
    l: l
  };
}

function hslShift(h, s, l) {
  return z => {
    let hsl = rgbToHsl(z.red, z.green, z.blue);
    return {
      ...z,
      ...hslToRgb((hsl.h + h + 1) % 1, hsl.s + s, hsl.l + l)
    }
  }
}

function setHue(hue) {
  return z => {
    let hsl = rgbToHsl(z.red, z.green, z.blue);
    return {
      ...z,
      ...hslToRgb(hue, hsl.s, hsl.l)
    }
  }
}

function setAlpha(alpha) {
  return z => {
    return {
      ...z,
      alpha: alpha
    }
  }
}

function setSaturation(saturation) {
  return z => {
    let hsl = rgbToHsl(z.red, z.green, z.blue);
    return {
      ...z,
      ...hslToRgb(hsl.h, saturation, hsl.l)
    }
  }
}

function setLightness(lightness) {
  return z => {
    let hsl = rgbToHsl(z.red, z.green, z.blue);
    return {
      ...z,
      ...hslToRgb(hsl.h, hsl.s, lightness)
    }
  }
}

function lerpRGB(r, g, b, weight = 0.5) {
  weight = Math.max(0, Math.min(1, weight));
  return z => {
    return {
      ...z,
      ...{
        red: (z.red * (1 - weight) + r * weight),
        green: (z.green * (1 - weight) + g * weight),
        blue: (z.blue * (1 - weight) + b * weight)
      }
    }
  }
}

function lerpHSL(h, s, l, weight = 0.5) {
  weight = Math.max(0, Math.min(1, weight));
  return z => {
    let hsl = rgbToHsl(z.red, z.green, z.blue);
    return {
      ...z,
      ...hslToRgb(
        Math.abs((hsl.h % 1) - (h % 1)) > 0.5 ?
        (((hsl.h + (hsl.h % 1 > h % 1 ? 0 : 1)) * (1 - weight) + (h + (hsl.h % 1 > h % 1 ? 1 : 0)) * weight) + 1) % 1 :
        (hsl.h % 1) * (1 - weight) + (h % 1) * weight,
        hsl.s * (1 - weight) + s * weight,
        hsl.l * (1 - weight) + l * weight
      )
    }
  }
}

function normalizeRGB(color) {
  let magnitude = Math.max(color.red, color.green, color.blue);
  return {
    red: color.red / magnitude,
    green: color.green / magnitude,
    blue: color.blue / magnitude
  }
}

function normalizeColor() {
  return z => {
    let magnitude = Math.max(z.red, z.green, z.blue);
    return {
      ...z,
      red: z.red / magnitude,
      green: z.green / magnitude,
      blue: z.blue / magnitude
    }
  }
}

function brightenRGB(color, amount) {
  return {
    red: Math.min(1, color.red * amount),
    green: Math.min(1, color.green * amount),
    blue: Math.min(1, color.blue * amount)
  }
}

function brighten(amount) {
  return z => {
    return {
      ...z,
      alpha: z.alpha * amount,
    }
  }
}

function repeatingGradient(colors) {
  if(colors.length < 1) {
    throw "not enough colors";
  }
  if(colors.length === 1) {
    let col = colors[0];
    if(colors[0].hasOwnProperty('h')) {
      col = hslToRgb(colors[0].h, colors[0].s, colors[0].l);
    }
    return z => {
      return {
        ...z,
        ...col
      }
    }
  }
  return z => {
    let at = (z.im - Math.floor(z.im)) * colors.length;
    return {
      ...z,
      ...lerp(colors[Math.floor(at) % colors.length], colors[Math.ceil(at) % colors.length], at % 1)
    }
  }
}

function gradient(colors) {
  if(colors.length < 1) {
    throw "not enough colors";
  }
  if(colors.length === 1) {
    let col = colors[0];
    if(colors[0].hasOwnProperty('h')) {
      col = hslToRgb(colors[0].h, colors[0].s, colors[0].l);
    }
    return z => {
      return {
        ...z,
        ...col
      }
    }
  }
  return z => {
    let at = Math.max(0, Math.min(1, z.im + 0.5)) * (colors.length - 1);
    return {
      ...z,
      ...lerp(colors[Math.floor(at) % colors.length], colors[Math.ceil(at) % colors.length], at % 1)
    }
  }
}

function gamma(gamma) {
  return z => {
    return {
      ...z,
      red: Math.pow(z.red, 1 / gamma),
      green: Math.pow(z.green, 1 / gamma),
      blue: Math.pow(z.blue, 1 / gamma)
    }
  }
}

function colorAdd(a, b) {
  return {
    red: a.red + b.red,
    green: a.green + b.green,
    blue: a.blue + b.blue
  }
}

function colorMult(a, b) {
  return {
    red: a.red * b.red,
    green: a.green * b.green,
    blue: a.blue * b.blue
  }
}

function colorMultScalar(a, scalar) {
  return {
    red: a.red * scalar,
    green: a.green * scalar,
    blue: a.blue * scalar
  }
}

function colorCos(a) {
  return {
    red: Math.cos(a.red),
    green: Math.cos(a.green),
    blue: Math.cos(a.blue)
  }
}

function paletteMod(a, b, c, d) {
  return z => {
    let f = colorAdd(a, colorMult(b, colorCos(colorMultScalar(colorAdd(colorMult(c, z), d), 2 * Math.PI))))
    return {
      ...z,
      red: f.red,
      green: f.green,
      blue: f.blue,
    }
  }
}

function rainbowCirc(x, y, _z, start) {
  return z => {
    return {
      ...color(colorHSL(
        Math.abs(Math.log(start + Math.hypot(z.re - x, z.im - y, z.z - _z))) % 1,
        1, 0.5))(z)
    }
  }
}

function rainbowCircAdd(x, y, _z, start) {
  return z => {
    let oh = rgbToHsl(z.red, z.green, z.blue).h;
    return {
      ...color(colorHSL(
        Math.abs(Math.log((start + Math.hypot(z.re - x, z.im - y, z.z - _z)) * oh)) % 1,
        1, 0.5))(z)
    }
  }
}

/* matrix math */
const IDENTITY_MATRIX = [
  1, 0, 0, 0,
  0, 1, 0, 0,
  0, 0, 1, 0,
  0, 0, 0, 1
];

function toRadian(a) {
  return a * DEGREE;
}

function fix3DZero(n) {
  return n === 0 ? Number.MAX_SAFE_INTEGER : n;
}

function applyMatrix(point, matrix) {
  let x = point[0],
    y = point[1],
    z = point[2];
  let hom = x * matrix[12] + y * matrix[13] + z * matrix[14] + matrix[15];
  return ([
    (x * matrix[0] + y * matrix[1] + z * matrix[2] + matrix[3]) / hom,
    (x * matrix[4] + y * matrix[5] + z * matrix[6] + matrix[7]) / hom,
    (x * matrix[8] + y * matrix[9] + z * matrix[10] + matrix[11]) / hom,
    1
  ]);
}

function multiplyMatrixAndPoint(matrix, point) {
  // Give a simple variable name to each part of the matrix, a column and row number
  let c0r0 = matrix[0],
    c1r0 = matrix[1],
    c2r0 = matrix[2],
    c3r0 = matrix[3];
  let c0r1 = matrix[4],
    c1r1 = matrix[5],
    c2r1 = matrix[6],
    c3r1 = matrix[7];
  let c0r2 = matrix[8],
    c1r2 = matrix[9],
    c2r2 = matrix[10],
    c3r2 = matrix[11];
  let c0r3 = matrix[12],
    c1r3 = matrix[13],
    c2r3 = matrix[14],
    c3r3 = matrix[15];

  // Now set some simple names for the point
  let x = point[0];
  let y = point[1];
  let z = point[2];
  let w = point[3];

  // Multiply the point against each part of the 1st column, then add together
  let resultX = (x * c0r0) + (y * c0r1) + (z * c0r2) + (w * c0r3);

  // Multiply the point against each part of the 2nd column, then add together
  let resultY = (x * c1r0) + (y * c1r1) + (z * c1r2) + (w * c1r3);

  // Multiply the point against each part of the 3rd column, then add together
  let resultZ = (x * c2r0) + (y * c2r1) + (z * c2r2) + (w * c2r3);

  // Multiply the point against each part of the 4th column, then add together
  let resultW = (x * c3r0) + (y * c3r1) + (z * c3r2) + (w * c3r3);

  return [resultX, resultY, resultZ, resultW];
}

function matrixMultiply(a, b) {
  let a00 = a[0],
    a01 = a[1],
    a02 = a[2],
    a03 = a[3];
  let a10 = a[4],
    a11 = a[5],
    a12 = a[6],
    a13 = a[7];
  let a20 = a[8],
    a21 = a[9],
    a22 = a[10],
    a23 = a[11];
  let a30 = a[12],
    a31 = a[13],
    a32 = a[14],
    a33 = a[15];

  // Cache only the current line of the second matrix
  let b0 = b[0],
    b1 = b[1],
    b2 = b[2],
    b3 = b[3];

  let out = new Array(16);
  out[0] = b0 * a00 + b1 * a10 + b2 * a20 + b3 * a30;
  out[1] = b0 * a01 + b1 * a11 + b2 * a21 + b3 * a31;
  out[2] = b0 * a02 + b1 * a12 + b2 * a22 + b3 * a32;
  out[3] = b0 * a03 + b1 * a13 + b2 * a23 + b3 * a33;

  b0 = b[4];
  b1 = b[5];
  b2 = b[6];
  b3 = b[7];
  out[4] = b0 * a00 + b1 * a10 + b2 * a20 + b3 * a30;
  out[5] = b0 * a01 + b1 * a11 + b2 * a21 + b3 * a31;
  out[6] = b0 * a02 + b1 * a12 + b2 * a22 + b3 * a32;
  out[7] = b0 * a03 + b1 * a13 + b2 * a23 + b3 * a33;

  b0 = b[8];
  b1 = b[9];
  b2 = b[10];
  b3 = b[11];
  out[8] = b0 * a00 + b1 * a10 + b2 * a20 + b3 * a30;
  out[9] = b0 * a01 + b1 * a11 + b2 * a21 + b3 * a31;
  out[10] = b0 * a02 + b1 * a12 + b2 * a22 + b3 * a32;
  out[11] = b0 * a03 + b1 * a13 + b2 * a23 + b3 * a33;

  b0 = b[12];
  b1 = b[13];
  b2 = b[14];
  b3 = b[15];
  out[12] = b0 * a00 + b1 * a10 + b2 * a20 + b3 * a30;
  out[13] = b0 * a01 + b1 * a11 + b2 * a21 + b3 * a31;
  out[14] = b0 * a02 + b1 * a12 + b2 * a22 + b3 * a32;
  out[15] = b0 * a03 + b1 * a13 + b2 * a23 + b3 * a33;
  return out;
}

function crossProduct(p1, p2) {
  return ([
    p1[1] * p2[2] - p1[2] * p2[1],
    p1[2] * p2[0] - p1[0] * p2[2],
    p1[0] * p2[1] - p1[1] * p2[0]
  ]);
}

function dotProduct(a, b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
}

function vectorTimes(m, n) {
  return [m[0] * n, m[1] * n, m[2] * n];
}

function vectorSum(m1, m2) {
  return [m1[0] + m2[0], m1[1] + m2[1], m1[2] + m2[2]];
}

function lerp4(a, b, t) {
  let ax = a[0];
  let ay = a[1];
  let az = a[2];
  let aw = a[3];
  let out = new Array(4);
  out[0] = ax + t * (b[0] - ax);
  out[1] = ay + t * (b[1] - ay);
  out[2] = az + t * (b[2] - az);
  out[3] = aw + t * (b[3] - aw);
  return out;
}

function normalize3(a) {
  let x = a[0];
  let y = a[1];
  let z = a[2];
  let len = x * x + y * y + z * z;
  if(len > 0) {
    len = 1 / Math.sqrt(len);
  }
  let out = new Array(3);
  out[0] = x * len;
  out[1] = y * len;
  out[2] = z * len;
  return out;
}

function normalOf3Points(p1, p2, p3) {
  let dir =
    crossProduct(
      [
        p1[0] - p2[0],
        p1[1] - p2[1],
        fix3DZero(p1[2]) - fix3DZero(p2[2])
      ],
      [
        p1[0] - p3[0],
        p1[1] - p3[1],
        fix3DZero(p1[2]) - fix3DZero(p3[2])
      ]
    );
  return normalize3(dir);
}

function matrixTranslate(x, y, z) {
  return ([
    1, 0, 0, x,
    0, 1, 0, y,
    0, 0, 1, z,
    0, 0, 0, 1,
  ]);
}

function matrixScale(x, y, z) {
  return ([
    x, 0, 0, 0,
    0, y, 0, 0,
    0, 0, z, 0,
    0, 0, 0, 1,
  ]);
}

function matrixRotateX(t) {
  return ([
    1, 0, 0, 0,
    0, Math.cos(t), -Math.sin(t), 0,
    0, Math.sin(t), Math.cos(t), 0,
    0, 0, 0, 1,
  ]);
}

function matrixRotateY(t) {
  return ([
    Math.cos(t), 0, Math.sin(t), 0,
    0, 1, 0, 0,
    -Math.sin(t), 0, Math.cos(t), 0,
    0, 0, 0, 1,
  ]);
}

function matrixRotateZ(t) {
  return ([
    Math.cos(t), -Math.sin(t), 0, 0,
    Math.sin(t), Math.cos(t), 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1,
  ]);
}

function matrixPerspectiveProjection(n, f) {
  return ([
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, (n + f) / n, -f,
    0, 0, 1 / n, 0,
  ]);
}

/* Quaternions */

function qAdd(a, b) {
  return {
    qx: a.qx + b.qx,
    qy: a.qy + b.qy,
    qz: a.qz + b.qz,
    qw: a.qw + b.qw
  };
}

function qAddScalar(a, scalar) {
  let s = {
    qx: scalar,
    qy: 0,
    qz: 0,
    qw: 0
  };
  return qAdd(a, s);
}

function qSub(a, b) {
  return {
    qx: a.qx - b.qx,
    qy: a.qy - b.qy,
    qz: a.qz - b.qz,
    qw: a.qw - b.qw
  };
}

function qMult(a, b) {
  let new0 = a.qx * b.qx - a.qy * b.qy - a.qz * b.qz - a.qw * b.qw;
  let new1 = a.qx * b.qy + a.qy * b.qx + a.qz * b.qw - a.qw * b.qz;
  let new2 = a.qx * b.qz - a.qy * b.qw + a.qz * b.qx + a.qw * b.qy;
  let new3 = a.qx * b.qw + a.qy * b.qz - a.qz * b.qy + a.qw * b.qx;
  return {
    qx: new0,
    qy: new1,
    qz: new2,
    qw: new3
  };
}

function qMultScalar(a, scalar) {
  let s = {
    qx: scalar,
    qy: 0,
    qz: 0,
    qw: 0
  };
  return qMult(a, s);
}

function qRecip(a) {
  let norm = (a.qx ** 2 + a.qy ** 2 + a.qz ** 2 + a.qw ** 2);
  return {
    qx: a.qx / norm,
    qy: -a.qy / norm,
    qz: -a.qz / norm,
    qw: -a.qw / norm
  };
}

function qDiv1(a, b) {
  return qMult(qRecip(b), a);
}

function qDiv2(a, b) {
  return qMult(a, qRecip(b));
}

function qMatNorm(mat) {
  let a = mat.q0;
  let b = mat.q1;
  let c = mat.q2;
  let d = mat.q3;

  let det = qSub(qMult(a, d), qMult(c, b));
  let denom = qMult(det, det);

  return {
    q0: qDiv1(a, denom),
    q1: qDiv1(b, denom),
    q2: qDiv1(c, denom),
    q3: qDiv1(d, denom)
  };
}

function qMatInv(mat) {
  let rm1 = {
    qx: -1,
    qy: 0,
    qz: 0,
    qw: 0
  };
  return {
    q0: mat.q3,
    q1: qMult(mat.q1, rm1),
    q2: qMult(mat.q2, rm1),
    q3: mat.q0
  };
}

function qConj(a) {
  return {
    qx: a.qx,
    qy: -a.qy,
    qz: -a.qz,
    qw: -a.qw
  };
}

function qNeg(a) {
  return {
    qx: -a.qx,
    qy: -a.qy,
    qz: -a.qz,
    qw: -a.qw
  };
}

function qExpZ(z) {
  let e = {
    qx: Math.exp(z.re),
    qy: 0,
    qz: 0,
    qw: 0
  };
  let v = {
    qx: 0,
    qy: z.im,
    qz: z.z,
    qw: 0
  };
  let absv = Math.hypot(z.im, z.z);

  let cv = {
    qx: Math.cos(absv),
    qy: 0,
    qz: 0,
    qw: 0
  };
  let sv = {
    qx: Math.sin(absv),
    qy: 0,
    qz: 0,
    qw: 0
  };
  let av = {
    qx: absv,
    qy: 0,
    qz: 0,
    qw: 0
  };

  return qMult(e, qAdd(cv, qMult(qDiv1(v, av), sv)));
}

function qExp(a) {
  let e = {
    qx: Math.exp(a.qx),
    qy: 0,
    qz: 0,
    qw: 0
  };
  let v = {
    qx: 0,
    qy: a.qy,
    qz: a.qz,
    qw: 0
  };
  let absv = Math.hypot(a.qy, a.qz);

  let cv = {
    qx: Math.cos(absv),
    qy: 0,
    qz: 0,
    qw: 0
  };
  let sv = {
    qx: Math.sin(absv),
    qy: 0,
    qz: 0,
    qw: 0
  };
  let av = {
    qx: absv,
    qy: 0,
    qz: 0,
    qw: 0
  };

  return qMult(e, qAdd(cv, qMult(qDiv1(v, av), sv)));
}

function qLogZ(z) {
  let absq = Math.hypot(z.re, z.im, z.z, 0);
  let v = {
    qx: 0,
    qy: z.im,
    qz: z.z,
    qw: 0
  };
  let absv = Math.hypot(z.im, z.z);

  let qln = {
    qx: Math.log(absq),
    qy: 0,
    qz: 0,
    qw: 0
  };
  let qac = {
    qx: Math.acos(z.re / absq),
    qy: 0,
    qz: 0,
    qw: 0
  };
  let qav = {
    qx: absv,
    qy: 0,
    qz: 0,
    qw: 0
  };

  return qAdd(qln, qMult(qDiv1(v, qav), qac))
}

function qLog(a) {
  let absq = Math.hypot(a.qx, a.qy, a.qz, 0);
  let v = {
    qx: 0,
    qy: a.qy,
    qz: a.qz,
    qw: 0
  };
  let absv = Math.hypot(a.qy, a.qz);

  let qln = {
    qx: Math.log(absq),
    qy: 0,
    qz: 0,
    qw: 0
  };
  let qac = {
    qx: Math.acos(a.qx / absq),
    qy: 0,
    qz: 0,
    qw: 0
  };
  let qav = {
    qx: absv,
    qy: 0,
    qz: 0,
    qw: 0
  };

  return qAdd(qln, qMult(qDiv1(v, qav), qac))
}

function qSinh(a) {
  let tmp = {
    qx: 2,
    qy: 0,
    qz: 0,
    qw: 0
  };
  return qSub((qDiv1(qExp(a), tmp)), (qDiv1(qExp(qNeg(a)), tmp)));
}

function qCosh(a) {
  let tmp = {
    qx: 2,
    qy: 0,
    qz: 0,
    qw: 0
  };
  return qAdd((qDiv1(qExp(qNeg(a)), tmp)), (qDiv1(qExp(a), tmp)));
}

function qTanh(a) {
  return qDiv1(qSinh(a), qCosh(a));
}

function qPowC(a, power) {
  let pow = {
    qx: power,
    qy: 0,
    qz: 0,
    qw: 0
  };
  return qExp(qMult(qLog(a), pow));
}

/* 3D */
function mobius3D(ar, ai, aj, ak, br, bi, bj, bk, cr, ci, cj, ck, dr, di, dj, dk, norm) {
  let preNorm0 = {
    qx: ar,
    qy: ai,
    qz: aj,
    qw: ak
  };
  let preNorm1 = {
    qx: br,
    qy: bi,
    qz: bj,
    qw: bk
  };
  let preNorm2 = {
    qx: cr,
    qy: ci,
    qz: cj,
    qw: ck
  };
  let preNorm3 = {
    qx: dr,
    qy: di,
    qz: dj,
    qw: dk
  };
  let preNorm = {
    q0: preNorm0,
    q1: preNorm1,
    q2: preNorm2,
    q3: preNorm3
  };
  let mobius = norm == 0 ? qMatNorm(preNorm) : preNorm;

  return z => {
    mobius = ak != 0 || bk != 0 || ck != 0 || dk != 0 ? (Math.random() > 0.5 ? mobius : qMatInv(mobius)) : mobius;
    let zin = {
      qx: z.re,
      qy: z.im,
      qz: z.z,
      qw: 0
    };
    let zout = qDiv1(qAdd(qMult(zin, mobius.q0), mobius.q1), qAdd(qMult(zin, mobius.q2), mobius.q3));

    return {
      ...z,
      re: zout.qx,
      im: zout.qy,
      z: zout.qz
    }
  }
}

function hypershift3D(x, y, _z) {
  let p = {
    qx: x,
    qy: y,
    qz: _z,
    qw: 0
  };
  let qr = {
    qx: 1,
    qy: 0,
    qz: 0,
    qw: 0
  };
  return z => {
    let q = {
      qx: z.re,
      qy: z.im,
      qz: z.z,
      qw: 0
    };
    let hs = qDiv1(qAdd(q, p), qAdd(qMult(q, qConj(p)), qr));
    return {
      ...z,
      re: hs.qx,
      im: hs.qy,
      z: hs.qz
    }
  }
}


function bubble3D() {
  return z => {
    let r = 4 / ((z.re ** 2 + z.im ** 2 + z.z ** 2) + 4)
    return {
      ...z,
      re: z.re * r,
      im: z.im * r,
      z: z.z * r
    }
  }
}

function julian3D(power) {
  let abspow = Math.abs(power);
  let powc = (1 / power - 1) * 0.5;
  return z => {
    let zz = z.z / abspow;
    let r2d = dot(z, z);
    let r = (r2d + zz * zz) ** powc;

    let r2 = r * r2d ** 0.5;
    let rnd = Math.floor(Math.random() * abspow);
    let ang = (Math.atan2(z.im, z.re) + 2 * Math.PI * rnd) / power;
    return {
      ...z,
      re: r2 * Math.cos(ang),
      im: r2 * Math.sin(ang),
      z: r * zz
    }
  }
}

function juliaq3D(power, divisor) {
  let invpow = divisor / power;
  let absinvpow = Math.abs(invpow);
  let halfinvpow = 0.5 * invpow - 0.5;
  let invpow2pi = 2 * Math.PI / power;
  return z => {
    let ang = Math.atan2(z.im, z.re) * invpow + Math.floor(Math.random() * 65536) * invpow2pi;

    let zz = z.z * absinvpow;
    let r2d = dot(z, z);
    let r = (r2d + ssqr(zz)) ** halfinvpow;
    let r2 = r * r2d ** 0.5;
    return {
      ...z,
      re: r2 * Math.cos(ang),
      im: r2 * Math.sin(ang),
      z: r * zz
    }
  }
}

function parabola3D() {
  return z => {
    let qz = {
      qx: z.re,
      qy: z.im,
      qz: z.z,
      qw: 0
    };
    return {
      ...z,
      ...qMult(qz, qz)
    }
  }
}

function sphereInv() {
  return z => {
    let r = 1 / (z.re ** 2 + z.im ** 2 + z.z ** 2);
    return {
      ...z,
      re: z.re * r,
      im: -z.im * r,
      z: -z.z * r
    }
  }
}

function trigCosh3D() {
  return z => {
    let qz = qCosh({
      qx: z.re,
      qy: z.im,
      qz: z.z,
      qw: 0
    });
    return {
      ...z,
      re: qz.qx,
      im: qz.qy,
      z: qz.qz
    }
  }
}

function trigExp3D() {
  return z => {
    let qe = qExpZ(z);
    return {
      ...z,
      re: qe.qx,
      im: qe.qy,
      z: qe.qz,
    }
  }
}

function trigLog3D() {
  return z => {
    let ql = qLogZ(z);
    return {
      ...z,
      re: ql.qx,
      im: ql.qy,
      z: ql.qz,
    }
  }
}

function trigSinh3D() {
  return z => {
    let qz = qSinh({
      qx: z.re,
      qy: z.im,
      qz: z.z,
      qw: 0
    });
    return {
      ...z,
      re: qz.qx,
      im: qz.qy,
      z: qz.qz
    }
  }
}

function trigTanh3D() {
  return z => {
    let qz = qTanh({
      qx: z.re,
      qy: z.im,
      qz: z.z,
      qw: 0
    });
    return {
      ...z,
      re: qz.qx,
      im: qz.qy,
      z: qz.qz
    }
  }
}

function unbubble3D() {
  return z => {
    let r = (z.re ** 2 + z.im ** 2 + z.z ** 2);
    let b = (Math.SQRT2 - Math.sqrt(2 - r)) / r;
    return {
      ...z,
      re: z.re * b,
      im: z.im * b,
      z: z.z * b
    }
  }
}

function matrix3D(matrix) {
  return z => {
    let result = applyMatrix([z.re, z.im, z.z], matrix);
    return {
      ...z,
      re: result[0],
      im: result[1],
      z: result[2],
    }
  }
}

function viewSphere(x, y, _z, radius) {
  let rst = reset();
  return z => {
    if(Math.hypot(z.re - x, z.im - y, z.z - _z) > radius) {
      return rst(z);
    }
    return z;
  }
}

function viewBox(x, y, _z, radius) {
  let rst = reset();
  return z => {
    if(z.re - x > radius || z.im - y > radius || z.z - _z > radius) {
      return rst(z);
    }
    return z;
  }
}

function getNormal(z) {
  if(z.z === 0 || z.re <= 0 || z.im <= 0 || z.re + 1 >= z.width || z.im + 1 >= z.height) {
    return [0, 0, 1];
  }
  let unit = 0.5 / Math.min(z.width, z.height);
  let n1 = normalOf3Points(
    [0, 0, z.z],
    [-unit, 0, z.zBuffer[(z.re - 1) + z.im * z.width]],
    [0, -unit, z.zBuffer[z.re + (z.im - 1) * z.width]]);
  let n2 = normalOf3Points(
    [0, 0, z.z],
    [0, -unit, z.zBuffer[z.re + (z.im - 1) * z.width]],
    [unit, 0, z.zBuffer[(z.re + 1) + z.im * z.width]]);
  let n3 = normalOf3Points(
    [0, 0, z.z],
    [0, unit, z.zBuffer[z.re + (z.im + 1) * z.width]],
    [-unit, 0, z.zBuffer[(z.re - 1) + z.im * z.width]]);
  let n4 = normalOf3Points(
    [0, 0, z.z],
    [unit, 0, z.zBuffer[(z.re + 1) + z.im * z.width]],
    [0, unit, z.zBuffer[z.re + (z.im + 1) * z.width]]);
  return normalize3([
    -(n1[0] + n2[0] + n3[0] + n4[0]),
    n1[1] + n2[1] + n3[1] + n4[1],
    n1[2] + n2[2] + n3[2] + n4[2]
  ]);
}

function blurSphere() {
  const [ax, ay] = getIrrationals(2);
  let u = Math.random();
  let v = Math.random();
  return z => {
    u = (u + ax) % 1;
    v = (v + ay) % 1;
    let theta = 2 * Math.PI * u;
    let phi = Math.acos(2 * v - 1);
    return {
      ...z,
      re: Math.sin(phi) * Math.cos(theta),
      im: Math.sin(phi) * Math.sin(theta),
      z: Math.cos(phi)
    }
  }
}

function blurCube() {
  const [ax, ay] = getIrrationals(2);
  let atx = Math.random();
  let aty = Math.random();
  return z => {
    atx = (atx + ax) % 1;
    aty = (aty + ay) % 1;
    let face = Math.random() * 6 >> 0;
    let x = atx - 0.5;
    let y = aty - 0.5;
    switch (face) {
      case 0:
        return {
          ...z,
          re: x,
            im: y,
            z: -0.5
        };
      case 1:
        return {
          ...z,
          re: x,
            im: y,
            z: 0.5
        };
      case 2:
        return {
          ...z,
          re: x,
            im: -0.5,
            z: y
        };
      case 3:
        return {
          ...z,
          re: x,
            im: 0.5,
            z: y
        };
      case 4:
        return {
          ...z,
          re: -0.5,
            im: x,
            z: y
        };
      case 5:
        return {
          ...z,
          re: 0.5,
            im: x,
            z: y
        };
    }
  }
}

function blurCubeVolume() {
  const [ax, ay, az] = getIrrationals(3);
  let atx = Math.random();
  let aty = Math.random();
  let atz = Math.random();
  return z => {
    atx = (atx + ax) % 1;
    aty = (aty + ay) % 1;
    atz = (atz + az) % 1;

    return {
      ...z,
      re: atx - 0.5,
      im: aty - 0.5,
      z: atz - 0.5
    };
  }
}

function blurSphereVolume() {
  const [ax, ay, az] = getIrrationals(3);
  let atx = Math.random();
  let aty = Math.random();
  let atz = Math.random();
  return z => {
    do {
      atx = (atx + ax) % 1;
      aty = (aty + ay) % 1;
      atz = (atz + az) % 1;
    } while((atx - 0.5) ** 2 + (aty - 0.5) ** 2 + (atz - 0.5) ** 2 > 0.25);

    return {
      ...z,
      re: atx * 2 - 1,
      im: aty * 2 - 1,
      z: atz * 2 - 1
    };
  }
}

function scale3D(s) {
  return z => {
    return {
      ...z,
      re: z.re * s,
      im: z.im * s,
      z: z.z * s
    }
  }
}

function scale3D3(x, y, _z) {
  return z => {
    return {
      ...z,
      re: z.re * x,
      im: z.im * y,
      z: z.z * _z
    }
  }
}

function translate3D(x, y, z) {
  return _z => {
    return {
      ..._z,
      re: _z.re + x,
      im: _z.im + y,
      z: _z.z + z
    }
  }
}

function rotate3D(pitch, roll, yaw) {
  let cosa = Math.cos(yaw * DEGREE);
  let sina = Math.sin(yaw * DEGREE);

  let cosb = Math.cos(pitch * DEGREE);
  let sinb = Math.sin(pitch * DEGREE);

  let cosc = Math.cos(roll * DEGREE);
  let sinc = Math.sin(roll * DEGREE);

  let Axx = cosa * cosb;
  let Axy = cosa * sinb * sinc - sina * cosc;
  let Axz = cosa * sinb * cosc + sina * sinc;

  let Ayx = sina * cosb;
  let Ayy = sina * sinb * sinc + cosa * cosc;
  let Ayz = sina * sinb * cosc - cosa * sinc;

  let Azx = -sinb;
  let Azy = cosb * sinc;
  let Azz = cosb * cosc;

  return z => {
    let px = z.re;
    let py = z.im;
    let pz = z.z;
    return {
      ...z,
      re: Axx * px + Axy * py + Axz * pz,
      im: Ayx * px + Ayy * py + Ayz * pz,
      z: Azx * px + Azy * py + Azz * pz,
    }
  }
}

function perspective3D() {
  let perspectiveMatrix = matrixMultiply(IDENTITY_MATRIX, matrixPerspectiveProjection(1, 0));
  return z => {
    let result = applyMatrix([z.re, z.im, z.z], perspectiveMatrix);
    if(z.z < 0) {
      return {
        ...z,
        re: Infinity,
        im: Infinity,
        z: 0
      }
    }
    return {
      ...z,
      re: result[0],
      im: result[1],
      z: Math.sqrt(Math.hypot(z.re, z.im, z.z))
    }
  }
}

function normalMap() {
  return z => {
    let normal = getNormal(z);
    return {
      ...z,
      red: (normal[0] + 1) / 2,
      green: (normal[1] + 1) / 2,
      blue: (normal[2] + 1) / 2,
    }
  }
}

function heightMap() {
  let low = Infinity;
  let high = -Infinity;
  return z => {
    if(z.z === 0 && z.red === 0 && z.green === 0 && z.blue === 0) {
      return {
        ...z,
        red: 0,
        blue: 0,
        green: 0
      };
    }
    if(z.z > high) {
      high = z.z;
    }
    if(z.z < low) {
      low = z.z;
    }
    let height = 1 - (z.z - low) / (high - low);
    return {
      ...z,
      red: height,
      green: height,
      blue: height,
    }
  }
}

function basicLighting(theta, diffuse) {
  let sinTheta = Math.sin(theta);
  let cosTheta = Math.cos(theta);

  return z => {
    let normal = getNormal(z);
    let brightness =
      Math.max(0, (normal[2] * sinTheta +
        normal[0] * cosTheta)) * (1 - diffuse) + diffuse;
    return {
      ...z,
      red: z.red * brightness,
      green: z.green * brightness,
      blue: z.blue * brightness,
    };
  };
}

function mist(startZ, halfLength, mistColor) {
  return z => {
    let brightness =
      z.z === 0 ? 0 :
      halfLength /
      (halfLength + Math.max(0, (z.z - startZ) / startZ));
    brightness = brightness * brightness;
    return {
      ...z,
      red: z.red * brightness + (1 - brightness) * mistColor.red,
      green: z.green * brightness + (1 - brightness) * mistColor.green,
      blue: z.blue * brightness + (1 - brightness) * mistColor.blue,
    }
  }
}

// TODO: maybe use textures instead
function environmentLight(x, y, z, lights) {
  let left = 1;
  let sum = 0;
  for(let i = 0; i < lights.length; i++) {
    let d = Math.hypot(x - lights[i][0][0], y + lights[i][0][1], z - lights[i][0][2]);
    if(d < lights[i][1] + Math.SQRT2) {
      let t = Math.min(1, (lights[i][1] + Math.SQRT2 - d) / Math.SQRT2);
      sum += left * t * lights[i][2];
      left -= left * t;
    }
  }
  return sum;
}

function oneLightEnvironment(x, y, z) {
  let lights = [
    [normalize3([0, 0, 1]), 0.2, 1],
  ];
  for(let i = 0; i < lights.length; i++) {
    if(Math.hypot(x - lights[i][0][0], y + lights[i][0][1], z - lights[i][0][2]) < lights[i][1]) {
      return {
        red: lights[i][2],
        green: lights[i][2],
        blue: lights[i][2]
      }
    }
  }
  return {
    red: 0.2,
    green: 0.2,
    blue: 0.2,
  }
}

function oneLightLights(x, y, z) {
  return environmentLight(x, y, z, [
    [normalize3([0, 0, 1]), 0.2, 0.9]
  ]);
}

function lightRoomEnvironment(x, y, z) {
  let lights = [
    [normalize3([3, -3, 1]), 0.4, 1],
    [normalize3([-3, -3, 1]), 0.4, 1],
  ];
  for(let i = 0; i < lights.length; i++) {
    if(Math.hypot(x - lights[i][0][0], y + lights[i][0][1], z - lights[i][0][2]) < lights[i][1]) {
      return {
        red: lights[i][2],
        green: lights[i][2],
        blue: lights[i][2]
      }
    }
  }
  let backdrop = normalize3([0, 1, -1]);
  if(y < -0.2 || Math.hypot(x - backdrop[0], y + backdrop[1], z - backdrop[2]) < 1.2) {
    let c = Math.max(0.2, y) * 0.5 + 0.3;
    return {
      red: c,
      green: c,
      blue: c
    };
  }
  if(y > 0.3) {
    return {
      red: -y * 0.1 + 0.2,
      green: -y * 0.1 + 0.2,
      blue: -y * 0.1 + 0.2
    }
  }
  let c = 0.2 + 0.1 - Math.abs(y + 0.05) * 0.4;
  return {
    red: c,
    green: c,
    blue: c
  }
}

function lightRoomLights(x, y, z) {
  return environmentLight(x, y, z, [
    [normalize3([0, 1, 0]), 0, 0.1],
    [normalize3([3, -3, 1]), 0.4, 1],
    [normalize3([-3, -3, 1]), 0.4, 1],
  ]);
}

function dayEnvironment(x, y, z) {
  let lights = [
    [normalize3([-1, -4, 1]), 0.02, 1],
    [normalize3([1.5, -2, -3]), 0.31, 0.9],
    [normalize3([3.75, -4.5, -0.5]), 0.56, 0.9],
    [normalize3([2, -1.75, 2]), 0.33, 0.9],
  ];
  for(let i = 0; i < lights.length; i++) {
    if(Math.hypot(x - lights[i][0][0], y + lights[i][0][1], z - lights[i][0][2]) < lights[i][1]) {
      return {
        red: lights[i][2],
        green: lights[i][2],
        blue: lights[i][2]
      }
    }
  }
  if(y > 0) {
    let c = (y - 0.2) * 0.2 + 0.8;
    return {
      red: c * 0.8,
      green: c * 0.9,
      blue: c * 1
    }
  }
  if(y > -0.1) {
    let c = (y - 0.2) * 0.2 + 0.7;
    return {
      red: c * 0.8,
      green: c * 0.8,
      blue: c * 1
    }
  }
  let c1 = -y;
  let c2 = 1 + y;
  return {
    red: c1 * 0.3 + c2 * 0.5,
    green: c1 * 0.7 + c2 * 0.4,
    blue: c1 * 0.2 + c2 * 0.7
  };
}

function dayLights(x, y, z) {
  return environmentLight(x, y, z, [
    [normalize3([-1, -4, 1]), 0.02, 1],
    [normalize3([0, 1, 0]), 0, 0.1],
    [normalize3([1.5, -2, -3]), 0.2, 0.4],
    [normalize3([3.75, -4.5, -0.5]), 0.3, 0.4],
    [normalize3([2, -1.75, 2]), 0.2, 0.4],
  ]);
}

function schlick(ior, normal) {
  let theta = dotProduct(
    [normal[0], normal[1], normal[2], 0],
    [0, 0, 1, 0]
  );
  let val = (1 - ior) / (1 + ior);
  let r0 = val * val;
  return (r0 + (1 - r0) * Math.pow(1 - theta, 5));
}

function inverseProjection(z, u, v) {
  let _z = z * z;
  return [
    _z / 1 * u,
    _z / 1 * v,
    _z
  ];
}

function _OLD_ambientOcclusion(s) {
  let weights = [];
  let maxP = 0;
  for(let i = -s; i <= s; i++) {
    weights[i + s] = [];
    for(let j = -s; j <= s; j++) {
      if(!(i === 0 && j === 0)) {
        let d = 1 / Math.sqrt(i * i + j * j);
        weights[i + s][j + s] = d;
        maxP += d;
      }
    }
  }
  return z => {
    let acc = 0;
    for(let i = -s; i <= s; i++) {
      for(let j = -s; j <= s; j++) {
        if(z.re + i < 0 || z.re + i >= z.width || z.im + j < 0 || z.im + j >= z.height) {
          acc += weights[i + s][j + s];
          continue;
        }
        if(!(i === 0 && j === 0)) {
          let d = weights[i + s][j + s];
          let sample = z.zBuffer[z.re + i + (z.im + j) * z.width];
          if(sample === 0 || sample > z.z) {
            acc += d;
          }
        }
      }
    }
    let c = acc / maxP;
    return {
      ...z,
      red: c,
      green: c,
      blue: c
    }
  }
}

function ambientOcclusion(steps, sz) {
  let slopeVectors = [
    [-1, 0],
    [1, 0],
    [0, -1],
    [0, 1],

    [-1, -1],
    [-1, 1],
    [1, -1],
    [1, 1],

    [-2, -1],
    [-1, -2],
    [-2, 1],
    [-1, 2],
    [2, -1],
    [1, -2],
    [2, 1],
    [1, 2],

    [-3, -2],
    [-3, -1],
    [-2, -3],
    [-1, -3],
    [-3, 2],
    [-3, 1],
    [-2, 3],
    [-1, 3],
    [3, -2],
    [3, -1],
    [2, -3],
    [1, -3],
    [3, 2],
    [3, 1],
    [2, 3],
    [1, 3],
    [-4, -3],
    [-4, -1],
    [-3, -4],
    [-1, -4],
    [-4, 3],
    [-4, 1],
    [-3, 4],
    [-1, 4],
    [4, -3],
    [4, -1],
    [3, -4],
    [1, -4],
    [4, 3],
    [4, 1],
    [3, 4],
    [1, 4],
    /*
    [-5,-4],
    [-5,-3],
    [-5,-2],
    [-5,-1],
    [-4,-5],
    [-3,-5],
    [-2,-5],
    [-1,-5],

    [-5,4],
    [-5,3],
    [-5,2],
    [-5,1],
    [-4,5],
    [-3,5],
    [-2,5],
    [-1,5],

    [5,-4],
    [5,-3],
    [5,-2],
    [5,-1],
    [4,-5],
    [3,-5],
    [2,-5],
    [1,-5],

    [5,4],
    [5,3],
    [5,2],
    [5,1],
    [4,5],
    [3,5],
    [2,5],
    [1,5],


    [-6,-5],
    [-6,-1],
    [-5,-6],
    [-1,-6],
    [-6,5],
    [-6,1],
    [-5,6],
    [-1,6],
    [6,-5],
    [6,-1],
    [5,-6],
    [1,-6],
    [6,5],
    [6,1],
    [5,6],
    [1,6],

    [-7,-6],
    [-7,-5],
    [-7,-4],
    [-7,-3],
    [-7,-2],
    [-7,-1],
    [-6,-7],
    [-5,-7],
    [-4,-7],
    [-3,-7],
    [-2,-7],
    [-1,-7],

    [-7,6],
    [-7,5],
    [-7,4],
    [-7,3],
    [-7,2],
    [-7,1],
    [-6,7],
    [-5,7],
    [-4,7],
    [-3,7],
    [-2,7],
    [-1,7],

    [7,-6],
    [7,-5],
    [7,-4],
    [7,-3],
    [7,-2],
    [7,-1],
    [6,-7],
    [5,-7],
    [4,-7],
    [3,-7],
    [2,-7],
    [1,-7],

    [7,6],
    [7,5],
    [7,4],
    [7,3],
    [7,2],
    [7,1],
    [6,7],
    [5,7],
    [4,7],
    [3,7],
    [2,7],
    [1,7],


    [10,1],
    [15,1],
    [17,1],
    [21,1],
    [11,10],
    [18,17],
    [1,10],
    [1,15],
    [1,17],
    [1,21],
    [10,11],
    [17,18],

    [-10,-1],
    [-15,-1],
    [-17,-1],
    [-21,-1],
    [-11,-10],
    [-18,-17],
    [-1,-10],
    [-1,-15],
    [-1,-17],
    [-1,-21],
    [-10,-11],
    [-17,-18],

    [10,-1],
    [15,-1],
    [17,-1],
    [21,-1],
    [11,-10],
    [18,-17],
    [1,-10],
    [1,-15],
    [1,-17],
    [1,-21],
    [10,-11],
    [17,-18],

    [-10,1],
    [-15,1],
    [-17,1],
    [-21,1],
    [-11,10],
    [-18,17],
    [-1,10],
    [-1,15],
    [-1,17],
    [-1,21],
    [-10,11],
    [-17,18],
    */
  ];
  for(let i = 0; i < slopeVectors.length; i++) {
    slopeVectors[i].push(Math.hypot(slopeVectors[i][0], slopeVectors[i][1]));
  }
  return z => {
    let s = sz / Math.min(z.width, z.height);
    if(z.z === 0 && z.red === 0 && z.green === 0 && z.blue === 0) {
      return {
        ...z,
        red: 0,
        blue: 0,
        green: 0
      };
    }

    let sum = 0;

    for(let j = 0; j < slopeVectors.length; j++) {
      let bestSlope = Number.MIN_SAFE_INTEGER;
      for(let i = 1; i < steps; i++) {
        let x = slopeVectors[j][0];
        let y = slopeVectors[j][1];
        let l = slopeVectors[j][2];
        let d = l * i * s;
        if(z.re + i * x < 0 || z.re + i * x >= z.width || z.im + i * y < 0 || z.im + i * y >= z.height) {
          i = steps;
          continue;
        }
        let N = z.re + i * x + (z.im + i * y) * z.width;
        if(z.zBuffer[N] === 0 && z.mainBuffer[0][N] === 0 && z.mainBuffer[1][N] === 0 && z.mainBuffer[2][N] === 0) {
          i = steps;
          continue;
        }
        let slope = (z.z - z.zBuffer[N]) / d;
        if(slope > bestSlope) {
          bestSlope = slope;
        }
      }
      sum += Math.atan(bestSlope);
    }

    sum /= Math.PI / 2 * slopeVectors.length;

    if(sum < 0) {
      sum = 0;
    }

    return {
      ...z,
      red: 1 - sum,
      blue: 1 - sum,
      green: 1 - sum
    };
  }
}

function _OLD_ambientOcclusionBig(s, step, depth) {
  let weights = [];
  let maxP = 0;
  for(let i = -s; i <= s; i++) {
    weights[i + s] = [];
    for(let j = -s; j <= s; j++) {
      if(!(i === 0 && j === 0)) {
        let d = 1 / Math.sqrt(i * i + j * j);
        weights[i + s][j + s] = d;
        maxP += d;
      }
    }
  }
  return z => {
    let acc = 0;
    for(let i = -s; i <= s; i++) {
      for(let j = -s; j <= s; j++) {
        if(z.re + i * step < 0 || z.re + i * step >= z.width || z.im + j * step < 0 || z.im + j * step >= z.height) {
          acc += weights[i + s][j + s];
          continue;
        }
        if(!(i === 0 && j === 0)) {
          let d = weights[i + s][j + s];
          let sample = z.zBuffer[z.re + i * step + (z.im + j * step) * z.width];
          if(sample === 0 || sample > z.z + depth) {
            acc += d;
          }
        }
      }
    }
    let c = acc / maxP;
    return {
      ...z,
      red: c,
      green: c,
      blue: c
    }
  }
}

function ambientOcclusionBig(steps, sz) {
  let slopeVectors = [
    [-1, 0],
    [1, 0],
    [0, -1],
    [0, 1],

    [-1, -1],
    [-1, 1],
    [1, -1],
    [1, 1],

    [-2, -1],
    [-1, -2],
    [-2, 1],
    [-1, 2],
    [2, -1],
    [1, -2],
    [2, 1],
    [1, 2],

    [-3, -2],
    [-3, -1],
    [-2, -3],
    [-1, -3],
    [-3, 2],
    [-3, 1],
    [-2, 3],
    [-1, 3],
    [3, -2],
    [3, -1],
    [2, -3],
    [1, -3],
    [3, 2],
    [3, 1],
    [2, 3],
    [1, 3],
    [-4, -3],
    [-4, -1],
    [-3, -4],
    [-1, -4],
    [-4, 3],
    [-4, 1],
    [-3, 4],
    [-1, 4],
    [4, -3],
    [4, -1],
    [3, -4],
    [1, -4],
    [4, 3],
    [4, 1],
    [3, 4],
    [1, 4],

    [-5, -4],
    [-5, -3],
    [-5, -2],
    [-5, -1],
    [-4, -5],
    [-3, -5],
    [-2, -5],
    [-1, -5],

    [-5, 4],
    [-5, 3],
    [-5, 2],
    [-5, 1],
    [-4, 5],
    [-3, 5],
    [-2, 5],
    [-1, 5],

    [5, -4],
    [5, -3],
    [5, -2],
    [5, -1],
    [4, -5],
    [3, -5],
    [2, -5],
    [1, -5],

    [5, 4],
    [5, 3],
    [5, 2],
    [5, 1],
    [4, 5],
    [3, 5],
    [2, 5],
    [1, 5],


    [-6, -5],
    [-6, -1],
    [-5, -6],
    [-1, -6],
    [-6, 5],
    [-6, 1],
    [-5, 6],
    [-1, 6],
    [6, -5],
    [6, -1],
    [5, -6],
    [1, -6],
    [6, 5],
    [6, 1],
    [5, 6],
    [1, 6],

    [-7, -6],
    [-7, -5],
    [-7, -4],
    [-7, -3],
    [-7, -2],
    [-7, -1],
    [-6, -7],
    [-5, -7],
    [-4, -7],
    [-3, -7],
    [-2, -7],
    [-1, -7],

    [-7, 6],
    [-7, 5],
    [-7, 4],
    [-7, 3],
    [-7, 2],
    [-7, 1],
    [-6, 7],
    [-5, 7],
    [-4, 7],
    [-3, 7],
    [-2, 7],
    [-1, 7],

    [7, -6],
    [7, -5],
    [7, -4],
    [7, -3],
    [7, -2],
    [7, -1],
    [6, -7],
    [5, -7],
    [4, -7],
    [3, -7],
    [2, -7],
    [1, -7],

    [7, 6],
    [7, 5],
    [7, 4],
    [7, 3],
    [7, 2],
    [7, 1],
    [6, 7],
    [5, 7],
    [4, 7],
    [3, 7],
    [2, 7],
    [1, 7],


    [10, 1],
    [15, 1],
    [17, 1],
    [21, 1],
    [11, 10],
    [18, 17],
    [1, 10],
    [1, 15],
    [1, 17],
    [1, 21],
    [10, 11],
    [17, 18],

    [-10, -1],
    [-15, -1],
    [-17, -1],
    [-21, -1],
    [-11, -10],
    [-18, -17],
    [-1, -10],
    [-1, -15],
    [-1, -17],
    [-1, -21],
    [-10, -11],
    [-17, -18],

    [10, -1],
    [15, -1],
    [17, -1],
    [21, -1],
    [11, -10],
    [18, -17],
    [1, -10],
    [1, -15],
    [1, -17],
    [1, -21],
    [10, -11],
    [17, -18],

    [-10, 1],
    [-15, 1],
    [-17, 1],
    [-21, 1],
    [-11, 10],
    [-18, 17],
    [-1, 10],
    [-1, 15],
    [-1, 17],
    [-1, 21],
    [-10, 11],
    [-17, 18],
  ];
  for(let i = 0; i < slopeVectors.length; i++) {
    slopeVectors[i].push(Math.hypot(slopeVectors[i][0], slopeVectors[i][1]));
  }
  return z => {
    let s = sz / Math.min(z.width, z.height);
    if(z.z === 0 && z.red === 0 && z.green === 0 && z.blue === 0) {
      return {
        ...z,
        red: 0,
        blue: 0,
        green: 0
      };
    }

    let sum = 0;

    for(let j = 0; j < slopeVectors.length; j++) {
      let bestSlope = Number.MIN_SAFE_INTEGER;
      for(let i = 1; i < steps; i++) {
        let x = slopeVectors[j][0];
        let y = slopeVectors[j][1];
        let l = slopeVectors[j][2];
        let d = l * i * s;
        if(z.re + i * x < 0 || z.re + i * x >= z.width || z.im + i * y < 0 || z.im + i * y >= z.height) {
          i = steps;
          continue;
        }
        let N = z.re + i * x + (z.im + i * y) * z.width;
        if(z.zBuffer[N] === 0 && z.mainBuffer[0][N] === 0 && z.mainBuffer[1][N] === 0 && z.mainBuffer[2][N] === 0) {
          i = steps;
          continue;
        }
        let slope = (z.z - z.zBuffer[N]) / d;
        if(slope > bestSlope) {
          bestSlope = slope;
        }
      }
      sum += Math.atan(bestSlope);
    }

    sum /= Math.PI / 2 * slopeVectors.length;

    if(sum < 0) {
      sum = 0;
    }

    return {
      ...z,
      red: 1 - sum,
      blue: 1 - sum,
      green: 1 - sum
    };
  }
}

function ambientOcclusion2(vectors, steps, sz, stepSize) {
  let slopeVectors = [];
  for(let i = 0; i < vectors; i++) {
    slopeVectors.push([Math.sin(i / vectors * Math.PI * 2), Math.cos(i / vectors * Math.PI * 2)]);
  }
  return z => {
    let s = sz / Math.min(z.width, z.height);
    if(z.z === 0 && z.red === 0 && z.green === 0 && z.blue === 0) {
      return {
        ...z,
        red: 0,
        blue: 0,
        green: 0
      };
    }

    let sum = 0;

    for(let j = 0; j < slopeVectors.length; j++) {
      let bestSlope = Number.MIN_SAFE_INTEGER;
      let bestD = steps * stepSize;
      let bestN = 1e10;
      for(let i = stepSize; i < steps * stepSize; i += stepSize) {
        let x = slopeVectors[j][0];
        let y = slopeVectors[j][1];

        let X = Math.round(z.re + i * x);
        let Y = Math.round(z.im + i * y);
        if(X < 0 || X >= z.width || Y < 0 || Y >= z.height) {
          i = Infinity;
          continue;
        }
        let N = X + Y * z.width;
        if(z.zBuffer[N] === 0 && z.mainBuffer[0][N] === 0 && z.mainBuffer[1][N] === 0 && z.mainBuffer[2][N] === 0) {
          i = Infinity;
          continue;
        }
        let d = Math.round(i * x) ** 2 + Math.round(i * y) ** 2;
        let slope = (z.z - z.zBuffer[N]) / d;
        if(slope > bestSlope) {
          bestSlope = slope;
          bestN = z.z - z.zBuffer[N];
          bestD = d;
        }
      }
      sum += Math.atan(bestN / (Math.sqrt(bestD) * s));
    }

    sum /= Math.PI / 2 * slopeVectors.length;

    if(sum < 0) {
      sum = 0;
    }

    return {
      ...z,
      red: 1 - sum,
      blue: 1 - sum,
      green: 1 - sum
    };
  }
}

function edgeDetection(sensitivityColor, sensitivityZ) {
  const neighbors = [
    [-1, -1],
    [-1, 0],
    [-1, 1],
    [0, -1],
    [0, 1],
    [1, -1],
    [1, 0],
    [1, 1]
  ];
  return z => {
    if(z.re === 0 || z.re === z.width - 1 || z.im === 0 || z.im === z.height - 1) {
      return {
        ...z,
        red: 1,
        green: 1,
        blue: 1
      };
    }

    let sum = 0;

    const at = z.re + z.im * z.width;
    const r = z.mainBuffer[0][at];
    const g = z.mainBuffer[1][at];
    const b = z.mainBuffer[2][at];
    const zz = z.zBuffer[at];

    for(let [x, y] of neighbors) {
      sum += Math.abs(z.mainBuffer[0][at + x + y * z.width] - r) * sensitivityColor;
      sum += Math.abs(z.mainBuffer[1][at + x + y * z.width] - g) * sensitivityColor;
      sum += Math.abs(z.mainBuffer[2][at + x + y * z.width] - b) * sensitivityColor;
      sum += Math.abs(z.zBuffer[at + x + y * z.width] - zz) * 40 * sensitivityZ;
    }

    sum /= 24;

    if(sum > 1) sum = 1;

    return {
      ...z,
      red: 1 - sum,
      blue: 1 - sum,
      green: 1 - sum
    };
  }
}

function reflect(vector, normal) {
  // v - 2 * (v dot n) * n
  return vectorSum(vector, vectorTimes(normal, -2 * (vector[0] * normal[0] + vector[1] * normal[1] + vector[2] * normal[2])));
}

function specularOrth(theta1, theta2, ior, environment, p) {
  let skyBoxLights = [lightRoomLights, dayLights, oneLightLights][environment];
  let rotation1 = rotate3D(theta1, 0, 0);
  let rotation2 = rotate3D(0, theta2, 0);

  return z => {
    if(z.z === 0 && z.red === 0 && z.green === 0 && z.blue === 0) {
      return {
        ...z,
        red: 0,
        blue: 0,
        green: 0
      };
    }
    let normal = getNormal(z);
    let reflected = reflect([0, 0, 1], [normal[0], normal[1], normal[2]]);
    let n = rotation1(rotation2({
      re: reflected[0],
      im: reflected[1],
      z: reflected[2]
    }));
    let brightness = Math.pow(skyBoxLights(-n.re, -n.im, -n.z), p);
    return {
      ...z,
      ...lerp({
          red: 0,
          green: 0,
          blue: 0
        }, {
          red: 1,
          green: 1,
          blue: 1
        },
        brightness * schlick(ior, normal))
    };
  }
}

function basicEnvironmentOrth(theta1, theta2, ior, environment) {
  let skyBox = [lightRoomEnvironment, dayEnvironment, oneLightEnvironment][environment];
  let skyBoxLights = [lightRoomLights, dayLights, oneLightLights][environment];
  let rotation1 = rotate3D(theta1, 0, 0);
  let rotation2 = rotate3D(0, theta2, 0);

  return z => {
    if(z.z === 0 && z.red === 0 && z.green === 0 && z.blue === 0) {
      return {
        ...z,
        red: 0,
        blue: 0,
        green: 0
      };
    }
    let normal = getNormal(z);
    let reflected = reflect([0, 0, 1], [normal[0], normal[1], normal[2]]);
    let norm = rotation1(rotation2({
      re: normal[0],
      im: normal[1],
      z: normal[2]
    }));
    let n = rotation1(rotation2({
      re: reflected[0],
      im: reflected[1],
      z: reflected[2]
    }));
    let brightness = skyBoxLights(norm.re, norm.im, norm.z);
    return {
      ...z,
      ...lerp(
        lerp({
            red: 0,
            green: 0,
            blue: 0
          },
          z, brightness),
        skyBox(-n.re, -n.im, -n.z),
        schlick(ior, normal)
      )
    };
  }
}

function advancedLightingOrth(theta1, theta2, ior, environment) {
  let p = 20;
  let specularWeight = 2;
  let skyBox = [lightRoomEnvironment, dayEnvironment, oneLightEnvironment][environment];
  let skyBoxLights = [lightRoomLights, dayLights, oneLightLights][environment];
  let rotation1 = rotate3D(theta1, 0, 0);
  let rotation2 = rotate3D(0, theta2, 0);

  return z => {
    let cameraSize = Math.min(z.width, z.height);
    let normal = getNormal(z);
    let reflected = reflect([0, 0, 1], [normal[0], normal[1], normal[2]]);
    let norm = rotation1(rotation2({
      re: normal[0],
      im: normal[1],
      z: normal[2]
    }));
    let n = rotation1(rotation2({
      re: reflected[0],
      im: reflected[1],
      z: reflected[2]
    }));
    let brightness = skyBoxLights(norm.re, norm.im, norm.z);
    if(normal[0] === 0 || normal[1] === 0 || !normal[2]) {
      let schlickValue = schlick(ior, normal);
      let specular = Math.pow(skyBoxLights(-n.re, -n.im, -n.z), p) * schlickValue * specularWeight;
      let reflections = lerp(
        lerp({
            red: 0,
            green: 0,
            blue: 0
          },
          z, brightness),
        skyBox(-n.re, -n.im, -n.z),
        schlickValue
      );
      return {
        ...z,
        red: reflections.red + specular,
        green: reflections.green + specular,
        blue: reflections.blue + specular,
      };
    }
    let amp = Math.sqrt(reflected[0] * reflected[0] + reflected[1] * reflected[1]);
    let step = [
      -reflected[0] / amp,
      reflected[1] / amp,
      reflected[2] / amp / cameraSize
    ];
    let a = Math.min(Math.abs(1 / step[0]), Math.abs(1 / step[1])) * 2;
    if(Math.abs(Math.abs(step[0]) - Math.SQRT1_2) < 0.1) {
      a *= 1.5;
    }
    let bounceNormalCheck = true;
    let _at = [z.re + step[0] * a, z.im + step[1] * a, z.z + step[2] * a];
    let steps = -1;
    while(true) {
      steps++;
      _at = [_at[0] + step[0], _at[1] + step[1], _at[2] + step[2]];
      let at = [Math.round(_at[0]), Math.round(_at[1]), _at[2]];
      if(at[0] < 0 || at[1] < 0 || at[0] + 1 >= z.width || at[1] + 1 >= z.height) {
        let schlickValue = schlick(ior, normal);
        let specular = Math.pow(skyBoxLights(-n.re, -n.im, -n.z), p) * schlickValue * specularWeight;
        let reflections = lerp(
          lerp({
              red: 0,
              green: 0,
              blue: 0
            },
            z, brightness),
          skyBox(-n.re, -n.im, -n.z),
          schlickValue
        );
        return {
          ...z,
          red: reflections.red + specular,
          green: reflections.green + specular,
          blue: reflections.blue + specular,
        };
      }
      let sample = z.zBuffer[at[0] + at[1] * z.width];

      let bounceNormal = getNormal({
        ...z,
        re: at[0],
        im: at[1],
        z: sample
      });
      if(
        sample !== 0 &&
        sample - at[2] > 0 === step[2] < 0 &&
        Math.abs(sample - at[2]) < Math.abs(step[2]) * 2
      ) {
        let newNorm = rotation1(rotation2({
          re: bounceNormal[0],
          im: bounceNormal[1],
          z: bounceNormal[2]
        }));
        let bounceReflected = reflect([-step[0], step[1], step[2]], [bounceNormal[0], bounceNormal[1], bounceNormal[2]]);
        let bounceN = rotation1(rotation2({
          re: bounceReflected[0],
          im: bounceReflected[1],
          z: bounceReflected[2]
        }));
        let bounceBrightness = skyBoxLights(newNorm.re, newNorm.im, newNorm.z);

        let bounceSchlickValue = schlick(ior, bounceNormal);
        let bounceSpecular = Math.pow(skyBoxLights(-bounceN.re, -bounceN.im, -bounceN.z), p) * bounceSchlickValue * specularWeight;

        let bounceReflections = lerp(
          lerp({
            red: 0,
            green: 0,
            blue: 0
          }, {
            red: z.mainBuffer[0][at[0] + at[1] * z.width],
            green: z.mainBuffer[1][at[0] + at[1] * z.width],
            blue: z.mainBuffer[2][at[0] + at[1] * z.width]
          }, bounceBrightness),
          skyBox(-bounceN.re, -bounceN.im, -bounceN.z),
          bounceSchlickValue
        );

        return {
          ...z,
          ...lerp(
            lerp({
                red: 0,
                green: 0,
                blue: 0
              },
              z, brightness), {
              red: bounceReflections.red + bounceSpecular,
              green: bounceReflections.green + bounceSpecular,
              blue: bounceReflections.blue + bounceSpecular,
            },
            schlick(ior, normal)
          )
        };
      }
    }
  }
}

function specular(theta1, theta2, ior, environment, p) {
  let skyBoxLights = [lightRoomLights, dayLights, oneLightLights][environment];
  let rotation1 = rotate3D(theta1, 0, 0);
  let rotation2 = rotate3D(0, theta2, 0);

  return z => {
    if(z.z === 0 && z.red === 0 && z.green === 0 && z.blue === 0) {
      return {
        ...z,
        red: 0,
        blue: 0,
        green: 0
      };
    }
    let normal = getNormal(z);
    let reflected = reflect([0, 0, -1], [normal[0], normal[1], -normal[2]]);
    let n = rotation1(rotation2({
      re: -reflected[0],
      im: -reflected[1],
      z: reflected[2]
    }));
    let brightness = Math.pow(skyBoxLights(n.re, n.im, n.z), p);
    return {
      ...z,
      ...lerp({
          red: 0,
          green: 0,
          blue: 0
        }, {
          red: 1,
          green: 1,
          blue: 1
        },
        brightness * schlick(ior, normal))
    };
  }
}

function basicEnvironment(theta1, theta2, ior, environment) {
  let skyBox = [lightRoomEnvironment, dayEnvironment, oneLightEnvironment][environment];
  let skyBoxLights = [lightRoomLights, dayLights, oneLightLights][environment];
  let rotation1 = rotate3D(theta1, 0, 0);
  let rotation2 = rotate3D(0, theta2, 0);

  return z => {
    if(z.z === 0 && z.red === 0 && z.green === 0 && z.blue === 0) {
      return {
        ...z,
        red: 0,
        blue: 0,
        green: 0
      };
    }
    let normal = getNormal(z);
    let n = rotation1(rotation2({
      re: normal[0],
      im: normal[1],
      z: normal[2]
    }));
    let brightness = skyBoxLights(n.re, n.im, n.z);
    return {
      ...z,
      ...lerp({
          red: 0,
          green: 0,
          blue: 0
        },
        lerp(
          z,
          skyBox(n.re, n.im, n.z),
          schlick(ior, normal)),
        brightness)
    };
  }
}

function advancedLighting(theta1, theta2, ior, environment) {
  let skyBox = [lightRoomEnvironment, dayEnvironment, oneLightEnvironment][environment];
  let skyBoxLights = [lightRoomLights, dayLights, oneLightLights][environment];
  let rotation1 = rotate3D(theta1, 0, 0);
  let rotation2 = rotate3D(0, theta2, 0);

  return z => {
    if(z.z === 0 && z.red === 0 && z.green === 0 && z.blue === 0) {
      return {
        ...z,
        red: 0,
        blue: 0,
        green: 0
      };
    }
    let cameraSize = Math.min(z.width, z.height);
    let normal = getNormal(z);
    let _normal = normal;
    let n = rotation1(rotation2({
      re: normal[0],
      im: normal[1],
      z: normal[2]
    }));
    let brightness = skyBoxLights(n.re, n.im, n.z);
    if(normal[0] === 0 || normal[1] === 0 || !normal[2]) {
      return {
        ...z,
        ...lerp({
            red: 0,
            green: 0,
            blue: 0
          },
          lerp(
            z,
            skyBox(n.re, n.im, n.z),
            schlick(ior, _normal)),
          brightness)
      };
    }
    let step = reflect([0, 0, 1], normal);
    let amp = Math.sqrt(step[0] * step[0] + step[1] * step[1]);
    step = [
      -step[0] / amp,
      -step[1] / amp,
      -step[2] / amp * // Math.log2(1+z.z*z.z)//z.z * z.z / amp
      (1 / z.z / cameraSize) //(-z.z*Math.abs(z.z / cameraSize / Math.SQRT2))
    ];
    //-z.z*Math.abs(z.z / amp / cameraSize / Math.SQRT2)]; // closeish
    //-amp / z.z / Math.min(z.width,z.height) * Math.SQRT2];
    let _at = [z.re, z.im, z.z];

    //if(z.re === 644 && (z.im === 344 || z.im === 223)){
    //  console.log(`${z.im}: ${z.z}`);
    //  console.log(z);
    //  console.log(amp);
    //  console.log(_normal);
    //}

    while(true) {
      _at = [_at[0] + step[0], _at[1] + step[1], _at[2] + step[2]];
      let at = [Math.round(_at[0]), Math.round(_at[1]), _at[2]];
      if(at[0] < 0 || at[1] < 0 || at[0] + 1 >= z.width || at[1] + 1 >= z.height) {
        return {
          ...z,
          ...lerp({
              red: 0,
              green: 0,
              blue: 0
            },
            lerp(
              z,
              skyBox(n.re, n.im, n.z),
              schlick(ior, _normal)),
            brightness)
        };
      }
      let sample = z.zBuffer[at[0] + at[1] * z.width];
      //if(sample !== 0 && sample > at[2] && Math.random() < 0.0001){
      //  console.log(`SAMPLE ${sample} > ${at[2]}\n${sample-at[2]} < ${normal[2]}`);
      //}
      if(sample !== 0 && sample > at[2] && Math.abs(sample - at[2]) < Math.abs(normal[2]) * 2) {
        bounceNormal = getNormal({
          ...z,
          re: at[0],
          im: at[1],
          z: sample
        });
        newN = rotation1(rotation2({
          re: bounceNormal[0],
          im: bounceNormal[1],
          z: bounceNormal[2]
        }));
        let bounceBrightness = skyBoxLights(n.re, n.im, n.z);
        return {
          ...z,
          ...lerp({
              red: 0,
              green: 0,
              blue: 0
            },
            lerp(
              z,
              lerp({
                  red: 0,
                  green: 0,
                  blue: 0
                },
                lerp({
                    red: z.mainBuffer[0][at[0] + at[1] * z.width],
                    green: z.mainBuffer[1][at[0] + at[1] * z.width],
                    blue: z.mainBuffer[2][at[0] + at[1] * z.width]
                  },
                  skyBox(newN.re, newN.im, newN.z),
                  schlick(ior, bounceNormal)),
                bounceBrightness),
              schlick(ior, _normal)),
            brightness)
        }
      }
    }
  }
}

/*images*/
async function getImage(URL) {
  let url = URL;
  if(url[0] === '/') {
    url = '/Net3arth/images' + url;
  }
  const response = await fetch(url);
  const blob = await response.blob();
  const bitmap = await createImageBitmap(blob);
  const offscreen = new OffscreenCanvas(bitmap.width, bitmap.height);
  const ctx = offscreen.getContext('2d');
  ctx.drawImage(bitmap, 0, 0, bitmap.width, bitmap.height);
  return ctx.getImageData(0, 0, bitmap.width, bitmap.height);
}

function setImage(url, x, y, w, h) {
  let img = false;
  let p = getImage(url);
  p.then(ans => {
    img = ans;
    promises.shift();
  });
  promises.push(p);
  return z => {
    if(!img) {
      return {
        ...z,
        red: 0,
        green: 0,
        blue: 0,
        alpha: 0,
        re: Infinity,
        im: Infinity
      };
    }
    let X = Math.floor((z.re - x) / w * img.width);
    let Y = Math.floor((z.im - y) / h * img.height);
    if(X < 0 || Y < 0 || X >= img.width || Y >= img.height) {
      return {
        ...z,
        red: 0,
        green: 0,
        blue: 0,
        alpha: 0,
      }
    }
    return {
      ...z,
      red: img.data[(X + Y * img.width) * 4] / 255,
      green: img.data[(X + Y * img.width) * 4 + 1] / 255,
      blue: img.data[(X + Y * img.width) * 4 + 2] / 255,
      alpha: img.data[(X + Y * img.width) * 4 + 3] / 255
    };
  }
}

function blurImage(url, x, y, w, h) {
  let img = false;
  let p = getImage(url);
  p.then(ans => {
    img = ans;
    promises.shift();
  });
  promises.push(p);
  const [ax, ay] = getIrrationals(2);
  let Xr = Math.random();
  let Yr = Math.random();
  return z => {
    if(!img) {
      return {
        ...z,
        red: 0,
        green: 0,
        blue: 0,
        alpha: 0,
        re: Infinity,
        im: Infinity
      };
    }
    let X, Y;
    do {
      Xr = (Xr + ax) % 1;
      Yr = (Yr + ay) % 1;
      X = Math.floor(Xr * img.width);
      Y = Math.floor(Yr * img.height);
    } while(img.data[(X + Y * img.width) * 4 + 3] === 0)
    return {
      ...z,
      re: x + Xr * w,
      im: y + Yr * h,
      red: img.data[(X + Y * img.width) * 4] / 255,
      green: img.data[(X + Y * img.width) * 4 + 1] / 255,
      blue: img.data[(X + Y * img.width) * 4 + 2] / 255,
      alpha: img.data[(X + Y * img.width) * 4 + 3] / 255
    };
  }
}

function matcap(url) {
  let img = false;
  let p = getImage(url);
  p.then(ans=>{
      img = ans;
      promises.shift();
    }
  );
  promises.push(p);
  return z=>{
    if(z.z === 0) {
      return z;
    }
    if (!img) {
      return {
        ...z,
        red: 0,
        green: 0,
        blue: 0,
        alpha: 0,
        re: Infinity,
        im: Infinity
      };
    }
    let normal = getNormal(z);
    let X = Math.floor((normal[0] + 1) / 2 * img.width);
    let Y = Math.floor((1 - normal[1]) / 2 * img.height);

    return {
      ...z,
      red: img.data[(X + Y * img.width) * 4] / 255,
      green: img.data[(X + Y * img.width) * 4 + 1] / 255,
      blue: img.data[(X + Y * img.width) * 4 + 2] / 255,
    };
  }
}

/*3D models*/
function decodeSTL(stl) {
  let file = Uint8Array.from(atob(stl), v => v.charCodeAt(0));

  let at = 84;

  let buf = new ArrayBuffer(4);
  let int8 = new Uint8Array(buf);
  let f32 = new Float32Array(buf);

  let n = file[80] + (file[81] << 8) + (file[82] << 16) + (file[83] << 24);

  let triangles = [];

  for(let i = 0; i < n; i++) {
    at += 12;

    int8[0] = file[at++];
    int8[1] = file[at++];
    int8[2] = file[at++];
    int8[3] = file[at++];
    let x1 = f32[0];
    int8[0] = file[at++];
    int8[1] = file[at++];
    int8[2] = file[at++];
    int8[3] = file[at++];
    let y1 = f32[0];
    int8[0] = file[at++];
    int8[1] = file[at++];
    int8[2] = file[at++];
    int8[3] = file[at++];
    let z1 = f32[0];
    int8[0] = file[at++];
    int8[1] = file[at++];
    int8[2] = file[at++];
    int8[3] = file[at++];

    let x2 = f32[0];
    int8[0] = file[at++];
    int8[1] = file[at++];
    int8[2] = file[at++];
    int8[3] = file[at++];
    let y2 = f32[0];
    int8[0] = file[at++];
    int8[1] = file[at++];
    int8[2] = file[at++];
    int8[3] = file[at++];
    let z2 = f32[0];
    int8[0] = file[at++];
    int8[1] = file[at++];
    int8[2] = file[at++];
    int8[3] = file[at++];

    let x3 = f32[0];
    int8[0] = file[at++];
    int8[1] = file[at++];
    int8[2] = file[at++];
    int8[3] = file[at++];
    let y3 = f32[0];
    int8[0] = file[at++];
    int8[1] = file[at++];
    int8[2] = file[at++];
    int8[3] = file[at++];
    let z3 = f32[0];

    at += 2;

    let triangle = [x1, y1, z1, x2, y2, z2, x3, y3, z3];

    triangles.push(triangle);
  }

  return triangles;
}

function stl(stl) {
  let triangles = decodeSTL(stl);
  let len = triangles.length;
  return z => {
    let t = Math.random() * len >> 0;
    return {
      ...z,
      ...triangle3D(...triangles[t])
    };
  }
}

function blurTetrahedron() {
  let v = [
    [1 / Math.sqrt(3), 1 / Math.sqrt(3), 1 / Math.sqrt(3)],
    [-1 / Math.sqrt(3), -1 / Math.sqrt(3), 1 / Math.sqrt(3)],
    [-1 / Math.sqrt(3), 1 / Math.sqrt(3), -1 / Math.sqrt(3)],
    [1 / Math.sqrt(3), -1 / Math.sqrt(3), -1 / Math.sqrt(3)],
  ];
  let triangles = [
    [...v[0], ...v[1], ...v[2]],
    [...v[0], ...v[1], ...v[3]],
    [...v[0], ...v[2], ...v[3]],
    [...v[1], ...v[2], ...v[3]],
  ];
  let len = triangles.length;
  return z => {
    return {
      ...z,
      ...triangle3D(...triangles[len * Math.random() >> 0])
    }
  }
}

function blurNormalCube() {
  let sqrt3 = 1 / Math.sqrt(3);
  return z => {
    let face = Math.random() * 6 >> 0;
    let x = Math.random() * sqrt3 * 2 - sqrt3;
    let y = Math.random() * sqrt3 * 2 - sqrt3;
    switch (face) {
      case 0:
        return {
          ...z,
          re: x,
            im: y,
            z: -sqrt3
        };
      case 1:
        return {
          ...z,
          re: x,
            im: y,
            z: sqrt3
        };
      case 2:
        return {
          ...z,
          re: x,
            im: -sqrt3,
            z: y
        };
      case 3:
        return {
          ...z,
          re: x,
            im: sqrt3,
            z: y
        };
      case 4:
        return {
          ...z,
          re: -sqrt3,
            im: x,
            z: y
        };
      case 5:
        return {
          ...z,
          re: sqrt3,
            im: x,
            z: y
        };
    }
  }
}

function blurOctahedron() {
  let v = [
    [1, 0, 0],
    [-1, 0, 0],
    [0, 1, 0],
    [0, -1, 0],
    [0, 0, 1],
    [0, 0, -1],
  ];
  let triangles = [
    [...v[0], ...v[2], ...v[4]],
    [...v[0], ...v[2], ...v[5]],
    [...v[0], ...v[3], ...v[4]],
    [...v[0], ...v[3], ...v[5]],
    [...v[1], ...v[2], ...v[4]],
    [...v[1], ...v[2], ...v[5]],
    [...v[1], ...v[3], ...v[4]],
    [...v[1], ...v[3], ...v[5]],
  ];
  let len = triangles.length;
  return z => {
    return {
      ...z,
      ...triangle3D(...triangles[len * Math.random() >> 0])
    }
  }
}

function blurIcosahedron() {
  let v = [
    [0, 0, -1],
    [0.7235999703407288, -0.5257200002670288, -0.4472149908542633],
    [-0.27638500928878784, -0.8506399989128113, -0.4472149908542633],
    [0.7235999703407288, 0.5257200002670288, -0.4472149908542633],
    [-0.8944249749183655, 0, -0.4472149908542633],
    [-0.27638500928878784, 0.8506399989128113, -0.4472149908542633],
    [0.8944249749183655, 0, 0.4472149908542633],
    [0.27638500928878784, -0.8506399989128113, 0.4472149908542633],
    [-0.7235999703407288, -0.5257200002670288, 0.4472149908542633],
    [-0.7235999703407288, 0.5257200002670288, 0.4472149908542633],
    [0.27638500928878784, 0.8506399989128113, 0.4472149908542633],
    [0, 0, 1]
  ];
  let triangles = [
    [...v[0], ...v[1], ...v[2]],
    [...v[0], ...v[1], ...v[3]],
    [...v[0], ...v[2], ...v[4]],
    [...v[0], ...v[4], ...v[5]],
    [...v[0], ...v[5], ...v[3]],
    [...v[1], ...v[3], ...v[6]],
    [...v[2], ...v[1], ...v[7]],
    [...v[4], ...v[2], ...v[8]],
    [...v[5], ...v[4], ...v[9]],
    [...v[3], ...v[5], ...v[10]],
    [...v[1], ...v[6], ...v[7]],
    [...v[2], ...v[7], ...v[8]],
    [...v[4], ...v[8], ...v[9]],
    [...v[5], ...v[9], ...v[10]],
    [...v[3], ...v[10], ...v[6]],
    [...v[7], ...v[6], ...v[11]],
    [...v[8], ...v[7], ...v[11]],
    [...v[9], ...v[8], ...v[11]],
    [...v[10], ...v[9], ...v[11]],
    [...v[6], ...v[10], ...v[11]]
  ];
  let len = triangles.length;
  return z => {
    return {
      ...z,
      ...triangle3D(...triangles[len * Math.random() >> 0])
    }
  }
}

function blurDodecahedron() {
  let v = [
    [-0.4910933282237111, -0.3567880810639798, -0.7946883717438183],
    [-0.49112925704843247, 0.35682693520010855, -0.7946660193643466],
    [0.1875920259326133, 0.5773398343994368, -0.794670069795526],
    [0.6070667364092035, 3.8720146025496536e-7, -0.7946644441966657],
    [0.1875959638518155, -0.5773508230692107, -0.794665869348377],
    [0.9822735907018106, 0.000012202298418618722, -0.18757158625675405],
    [0.7946518428552186, -0.5773298583374579, 0.1876724719963161],
    [0.30355513927219124, -0.9341778707812964, -0.18757751188755362],
    [-0.30355170765688644, -0.9341795959649469, 0.1875791245592269],
    [-0.7946531179909603, -0.5773277956178757, -0.18767024050876815],
    [-0.9822743407816587, 0.000011500861138844604, 0.18756961729715294],
    [-0.7946348910506528, 0.5773922649808152, -0.18759153838071205],
    [-0.30350175233900667, 0.9341631692162747, 0.18763485549193665],
    [0.30350152731505226, 0.9341629441923204, -0.1876356993317657],
    [0.794634816042668, 0.5773923774927924, 0.18759163214069305],
    [0.49109584099120207, -0.3568027076210167, 0.7946815460172011],
    [-0.1875955325559029, -0.5773510480931651, 0.7946657193324074],
    [-0.6070670739451351, 1.170227170155741e-7, 0.7946642941806961],
    [-0.18759144462073105, 0.5773415970870798, 0.7946690946917236],
    [0.49113120725603737, 0.35682723523204773, 0.7946643691886809]
  ];
  let triangles = [
    [...v[0], ...v[1], ...v[2]],
    [...v[2], ...v[3], ...v[4]],
    [...v[2], ...v[4], ...v[0]],
    [...v[4], ...v[3], ...v[5]],
    [...v[5], ...v[6], ...v[7]],
    [...v[5], ...v[7], ...v[4]],
    [...v[0], ...v[4], ...v[7]],
    [...v[7], ...v[8], ...v[9]],
    [...v[7], ...v[9], ...v[0]],
    [...v[1], ...v[0], ...v[9]],
    [...v[9], ...v[10], ...v[11]],
    [...v[9], ...v[11], ...v[1]],
    [...v[2], ...v[1], ...v[11]],
    [...v[11], ...v[12], ...v[13]],
    [...v[11], ...v[13], ...v[2]],
    [...v[14], ...v[5], ...v[3]],
    [...v[3], ...v[2], ...v[13]],
    [...v[3], ...v[13], ...v[14]],
    [...v[7], ...v[6], ...v[15]],
    [...v[15], ...v[16], ...v[8]],
    [...v[15], ...v[8], ...v[7]],
    [...v[9], ...v[8], ...v[16]],
    [...v[16], ...v[17], ...v[10]],
    [...v[16], ...v[10], ...v[9]],
    [...v[11], ...v[10], ...v[17]],
    [...v[17], ...v[18], ...v[12]],
    [...v[17], ...v[12], ...v[11]],
    [...v[13], ...v[12], ...v[18]],
    [...v[18], ...v[19], ...v[14]],
    [...v[18], ...v[14], ...v[13]],
    [...v[5], ...v[14], ...v[19]],
    [...v[19], ...v[15], ...v[6]],
    [...v[19], ...v[6], ...v[5]],
    [...v[17], ...v[16], ...v[15]],
    [...v[15], ...v[19], ...v[18]],
    [...v[15], ...v[18], ...v[17]]
  ];
  let len = triangles.length;
  return z => {
    return {
      ...z,
      ...triangle3D(...triangles[len * Math.random() >> 0])
    }
  }
}

function torus(r1, r2) {
  return z => {
    let t1 = Math.random() * Math.PI * 2;
    let t2 = Math.random() * Math.PI * 2;
    return {
      ...z,
      re: Math.cos(t1) * (r1 + Math.cos(t2) * r2),
      im: r2 * Math.sin(t2),
      z: Math.sin(t1) * (r1 + Math.cos(t2) * r2),
    };
  };
}

function blurTeapot() {
  return stl("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");
}

function umbilicTorus(mode) {
  return z => {
    let uv = mode == 0 ? C(Math.random()*2*Math.PI-Math.PI, Math.random()*2*Math.PI-Math.PI) : C(z.re * 2 * Math.PI, z.im * 2 * Math.PI);
    let q = (7 + Math.cos(uv.re/3-2*uv.im) + 2 * Math.cos(uv.re / 3 + uv.im));
    return {
      ...z,
      re: Math.sin(uv.re) * q * 0.25,
      im: Math.cos(uv.re) * q * 0.25,
      z: (Math.sin(uv.re/3 - 2*uv.im) + 2 * Math.sin(uv.re/3+uv.im)) * 0.25
    }
  }
}

/*dither*/
function bayerMatrixHelper(m) {
  var ans = [];
  var w = m.length;
  for(var i = 0; i < w; i++) {
    ans[i] = [];
    ans[i + w] = [];
    for(var j = 0; j < w; j++) {
      ans[i + 0][j + 0] = m[i][j] * 4;
      ans[i + w][j + w] = m[i][j] * 4 + 1;
      ans[i + 0][j + w] = m[i][j] * 4 + 2;
      ans[i + w][j + 0] = m[i][j] * 4 + 3;
    }
  }
  return ans;
}

function bayerMatrix(n) {
  var ans = [
    [0, 2],
    [3, 1]
  ];
  for(var i = 1; i < n; i++) {
    ans = bayerMatrixHelper(ans);
  }
  var l = ans.length;
  for(var i = 0; i < l; i++) {
    for(var j = 0; j < l; j++) {
      ans[i][j] = (ans[i][j] + 0.5) / (l * l);
    }
  }
  return ans;
}

function dither(s) {
  let ditherMatrix = bayerMatrix(s);
  let matrixSize = ditherMatrix.length;
  return z => {
    return {
      ...z,
      red: (z.red * 256 >> 0) / 256 + (ditherMatrix[z.re % matrixSize][z.im % matrixSize] > (z.red * 256 % 1) ? 0 : 1 / 256),
      green: (z.green * 256 >> 0) / 256 + (ditherMatrix[z.re % matrixSize][z.im % matrixSize] > (z.green * 256 % 1) ? 0 : 1 / 256),
      blue: (z.blue * 256 >> 0) / 256 + (ditherMatrix[z.re % matrixSize][z.im % matrixSize] > (z.blue * 256 % 1) ? 0 : 1 / 256),
    }
  }
}

/* descriptions*/
const BUILT_IN_TRANSFORMS = {
  identity: identity,
  draw: draw,
  dither: dither,
  //3D models
  stl: stl,
  blurTetrahedron: blurTetrahedron,
  blurNormalCube: blurNormalCube,
  blurOctahedron: blurOctahedron,
  blurIcosahedron: blurIcosahedron,
  blurDodecahedron: blurDodecahedron,
  blurTeapot: blurTeapot,
  torus: torus,
  umbilicTorus: umbilicTorus,
  //shaders
  shaderPass: shaderPass,
  normalizeColors: normalizeColors,
  rainbowCirc: rainbowCirc,
  rainbowCircAdd: rainbowCircAdd,
  paletteMod: paletteMod,
  gamma: gamma,
  ambientOcclusion: ambientOcclusion,
  ambientOcclusion2: ambientOcclusion2,
  ambientOcclusionBig: ambientOcclusionBig,
  edgeDetection: edgeDetection,
  specular: specular,
  specularOrth: specularOrth,
  normalMap: normalMap,
  heightMap: heightMap,
  matcap: matcap,
  basicLighting: basicLighting,
  basicEnvironment: basicEnvironment,
  basicEnvironmentOrth: basicEnvironmentOrth,
  advancedLighting: advancedLighting,
  advancedLightingOrth: advancedLightingOrth,
  mist: mist,
  //3D transforms
  mobius3D: mobius3D,
  hypershift3D: hypershift3D,
  bubble3D: bubble3D,
  julian3D: julian3D,
  juliaq3D: juliaq3D,
  parabola3D: parabola3D,
  sphereInv: sphereInv,
  trigCosh3D: trigCosh3D,
  trigExp3D: trigExp3D,
  trigLog3D: trigLog3D,
  trigSinh3D: trigSinh3D,
  trigTanh3D: trigTanh3D,
  unbubble3D: unbubble3D,
  matrix3D: matrix3D,
  blurSphere: blurSphere,
  blurSphereVolume: blurSphereVolume,
  blurCube: blurCube,
  blurCubeVolume: blurCubeVolume,
  scale3D: scale3D,
  scale3D3: scale3D3,
  translate3D: translate3D,
  rotate3D: rotate3D,
  perspective3D: perspective3D,
  viewBox: viewBox,
  viewSphere: viewSphere,
  //Images
  blurImage: blurImage,
  setImage: setImage,
  //2D transforms
  reset: reset,
  arcsinh: arcsinh,
  arctanh: arctanh,
  bent: bent,
  blurCircle: blurCircle,
  blurGasket: blurGasket,
  blurGaussian: blurGaussian,
  blurNgon: blurNgon,
  blurSine: blurSine,
  blurSquare: blurSquare,
  blurTriangle: blurTriangle,
  bTransform: bTransform,
  bubble: bubble,
  circleInv: circleInv,
  cpow: cpow,
  cylinder: cylinder,
  dc_poincareDisc: dc_poincareDisc,
  disc: disc,
  dragon: dragon,
  ePush: ePush,
  eRotate: eRotate,
  flipX: flipX,
  flipY: flipY,
  hypershape: hypershape,
  hypershift: hypershift,
  hypertile3: hypertile3,
  jac_cd: jac_cd,
  jac_cn: jac_cn,
  jac_dn: jac_dn,
  jac_elk: jac_elk,
  jac_sn: jac_sn,
  julian: julian,
  juliaq: juliaq,
  juliascope: juliascope,
  mobius: mobius,
  multiMobius: multiMobius,
  murl2: murl2,
  nSplit: nSplit,
  pdj: pdj,
  pointSymmetry: pointSymmetry,
  rotate: rotate,
  scale: scale,
  scale2: scale2,
  schwarzChristoffelMap: schwarzChristoffelMap,
  schwarzChristoffelInverseMap: schwarzChristoffelInverseMap,
  schwarzTriangle: schwarzTriangle,
  sinusoidal: sinusoidal,
  skew: skew,
  smartcrop: smartcrop,
  smartshape: smartshape,
  splits: splits,
  tileHelp: tileHelp,
  tileLog: tileLog,
  translate: translate,
  trigCosh: trigCosh,
  trigExp: trigExp,
  trigLog: trigLog,
  trigSinh: trigSinh,
  trigTanh: trigTanh,
  unbubble: unbubble,
  weierstrassElliptic: weierstrassElliptic,
  //color transfomrs
  brighten: brighten,
  color: color,
  gradient: gradient,
  repeatingGradient: repeatingGradient,
  hslShift: hslShift,
  lerpColor: lerpColor,
  lerpHSL: lerpHSL,
  lerpRGB: lerpRGB,
  normalizeColor: normalizeColor,
  setHue: setHue,
  setAlpha: setAlpha,
  setSaturation: setSaturation,
  setLightness: setLightness,
};

const BUILT_IN_TRANSFORMS_PARAMS = {
  identity: [],
  draw: [],
  dither: [{
    name: "matrix size",
    type: "number",
    default: 5
  }],
  //3D models
  stl: [{
    name: "stl",
    type: "string",
    default: ""
  }],
  blurTetrahedron: [],
  blurNormalCube: [],
  blurOctahedron: [],
  blurIcosahedron: [],
  blurDodecahedron: [],
  blurTeapot: [],
  torus: [{
    name: "Major Radius",
    type: "number",
    default: 0.7
  }, {
    name: "Minor Radius",
    type: "number",
    default: 0.3
  }],
  umbilicTorus: [{
    name: "mode",
    type: "number",
    default: 0
  }],
  //shaders
  shaderPass: [{
    name: "transform",
    type: "function",
    default: "heightMap()"
  }],
  normalizeColors: [],
  rainbowCirc: [{
      name: "x",
      type: "number",
      default: 0
    },
    {
      name: "y",
      type: "number",
      default: 0
    },
    {
      name: "z",
      type: "number",
      default: 0
    },
    {
      name: "start",
      type: "number",
      default: 1
    }
  ],
  rainbowCircAdd: [{
      name: "x",
      type: "number",
      default: 0
    },
    {
      name: "y",
      type: "number",
      default: 0
    },
    {
      name: "z",
      type: "number",
      default: 0
    },
    {
      name: "start",
      type: "number",
      default: 1
    }
  ],
  paletteMod: [{
      name: "a",
      type: "object",
      default: {
        red: 0.5,
        green: 0.5,
        blue: 0.5
      }
    },
    {
      name: "b",
      type: "object",
      default: {
        red: 0.5,
        green: 0.5,
        blue: 0.5
      }
    },
    {
      name: "c",
      type: "object",
      default: {
        red: 0.5,
        green: 0.5,
        blue: 0.5
      }
    },
    {
      name: "d",
      type: "object",
      default: {
        red: 0.5,
        green: 0.5,
        blue: 0.5
      }
    }
  ],
  gamma: [{
    name: "gamma",
    type: "number",
    default: 2.2
  }],
  ambientOcclusion: [{
    name: "steps",
    type: "number",
    default: 10
  }, {
    name: "size",
    type: "number",
    default: 1
  }],
  ambientOcclusion2: [{
    name: "vectors",
    type: "number",
    default: 48
  }, {
    name: "steps",
    type: "number",
    default: 10
  }, {
    name: "size",
    type: "number",
    default: 1
  }, {
    name: "stepSize",
    type: "number",
    default: 1
  }],
  ambientOcclusionBig: [{
    name: "steps",
    type: "number",
    default: 10
  }, {
    name: "size",
    type: "number",
    default: 1
  }],
  edgeDetection: [{
    name: "sensitivityColor",
    type: "number",
    default: 1
  }, {
    name: "sensitivityZ",
    type: "number",
    default: 1
  }],
  specular: [{
    name: "theta",
    type: "number",
    default: "20*DEGREE",
  }, {
    name: "theta",
    type: "number",
    default: "30*DEGREE",
  }, {
    name: "index of refraction",
    type: "number",
    default: 0.1,
  }, {
    name: "environment",
    type: "number",
    default: 0,
  }, {
    name: "power",
    type: "number",
    default: 20,
  }],
  specularOrth: [{
    name: "theta",
    type: "number",
    default: "20*DEGREE",
  }, {
    name: "theta",
    type: "number",
    default: "30*DEGREE",
  }, {
    name: "index of refraction",
    type: "number",
    default: 0.1,
  }, {
    name: "environment",
    type: "number",
    default: 0,
  }, {
    name: "power",
    type: "number",
    default: 20,
  }],
  normalMap: [],
  heightMap: [],
  matcap: [
    {
      name: "url",
      type: "string",
      default: "'/gold.png'"
    }
  ],
  basicLighting: [{
    name: "theta",
    type: "number",
    default: "30*DEGREE",
  }, {
    name: "diffuse",
    type: "number",
    default: 0.3,
  }],
  basicEnvironment: [{
    name: "theta",
    type: "number",
    default: "20*DEGREE",
  }, {
    name: "theta",
    type: "number",
    default: "30*DEGREE",
  }, {
    name: "index of refraction",
    type: "number",
    default: 1.4,
  }, {
    name: "environment",
    type: "number",
    default: 0,
  }],
  basicEnvironmentOrth: [{
    name: "theta",
    type: "number",
    default: "20*DEGREE",
  }, {
    name: "theta",
    type: "number",
    default: "30*DEGREE",
  }, {
    name: "index of refraction",
    type: "number",
    default: 1.4,
  }, {
    name: "environment",
    type: "number",
    default: 0,
  }],
  advancedLighting: [{
    name: "theta",
    type: "number",
    default: "20*DEGREE",
  }, {
    name: "theta",
    type: "number",
    default: "30*DEGREE",
  }, {
    name: "index of refraction",
    type: "number",
    default: 1.4,
  }, {
    name: "environment",
    type: "number",
    default: 0,
  }],
  advancedLightingOrth: [{
    name: "theta",
    type: "number",
    default: "20*DEGREE",
  }, {
    name: "theta",
    type: "number",
    default: "30*DEGREE",
  }, {
    name: "index of refraction",
    type: "number",
    default: 1.4,
  }, {
    name: "environment",
    type: "number",
    default: 0,
  }],
  mist: [{
    name: "startZ",
    type: "number",
    default: 1.5,
  }, {
    name: "halfLength",
    type: "number",
    default: 0.5,
  }, {
    name: "mistColor",
    type: "object",
    default: 'colorRGB(0.7, 0.7, 0.8)'
  }],
  //3D transforms
  mobius3D: [{
      name: "ar",
      type: "number",
      default: 1
    },
    {
      name: "ai",
      type: "number",
      default: 0
    },
    {
      name: "aj",
      type: "number",
      default: 0
    },
    {
      name: "ak",
      type: "number",
      default: 0
    },
    {
      name: "br",
      type: "number",
      default: 0
    },
    {
      name: "bi",
      type: "number",
      default: 0
    },
    {
      name: "bj",
      type: "number",
      default: 0
    },
    {
      name: "bk",
      type: "number",
      default: 0
    },
    {
      name: "cr",
      type: "number",
      default: 0
    },
    {
      name: "ci",
      type: "number",
      default: 0
    },
    {
      name: "cj",
      type: "number",
      default: 0
    },
    {
      name: "ck",
      type: "number",
      default: 0
    },
    {
      name: "dr",
      type: "number",
      default: 1
    },
    {
      name: "di",
      type: "number",
      default: 0
    },
    {
      name: "dj",
      type: "number",
      default: 0
    },
    {
      name: "dk",
      type: "number",
      default: 0
    },
    {
      name: "norm",
      type: "number",
      default: 1
    }
  ],
  hypershift3D: [{
      name: "x",
      type: "number",
      default: 0
    },
    {
      name: "y",
      type: "number",
      default: 0
    },
    {
      name: "z",
      type: "number",
      default: 0
    }
  ],
  matrix3D: [{
    name: "matrix",
    type: "array",
    default: IDENTITY_MATRIX
  }],
  bubble3D: [],
  julian3D: [{
    name: "power",
    type: "number",
    default: 1
  }],
  juliaq3D: [{
      name: "power",
      type: "number",
      default: 1
    },
    {
      name: "divisor",
      type: "number",
      default: 1
    }
  ],
  parabola3D: [],
  sphereInv: [],
  trigCosh3D: [],
  trigExp3D: [],
  trigLog3D: [],
  trigSinh3D: [],
  trigTanh3D: [],
  unbubble3D: [],
  blurSphere: [],
  blurSphereVolume: [],
  blurCube: [],
  blurCubeVolume: [],
  scale3D: [{
    name: "scale",
    type: "number",
    default: 1
  }],
  scale3D3: [{
    name: "x",
    type: "number",
    default: 1
  }, {
    name: "y",
    type: "number",
    default: 1
  }, {
    name: "z",
    type: "number",
    default: 1
  }],
  translate3D: [{
    name: "X",
    type: "number",
    default: 0
  }, {
    name: "Y",
    type: "number",
    default: 0
  }, {
    name: "Z",
    type: "number",
    default: 0
  }],
  rotate3D: [{
    name: "Pitch",
    type: "number",
    default: 0
  }, {
    name: "Roll",
    type: "number",
    default: 0
  }, {
    name: "Yaw",
    type: "number",
    default: 0
  }],
  perspective3D: [],
  viewBox: [{
    name: "x",
    type: "number",
    default: 0
  }, {
    name: "y",
    type: "number",
    default: 0
  }, {
    name: "z",
    type: "number",
    default: 0
  }, {
    name: "radius",
    type: "number",
    default: 1
  }, ],
  viewSphere: [{
    name: "x",
    type: "number",
    default: 0
  }, {
    name: "y",
    type: "number",
    default: 0
  }, {
    name: "z",
    type: "number",
    default: 0
  }, {
    name: "radius",
    type: "number",
    default: 1
  }, ],
  //Images
  blurImage: [{
    name: "url",
    type: "string",
    default: "'https://upload.wikimedia.org/wikipedia/en/thumb/8/80/Wikipedia-logo-v2.svg/1200px-Wikipedia-logo-v2.svg.png'"
  }, {
    name: "x",
    type: "number",
    default: -0.5
  }, {
    name: "y",
    type: "number",
    default: -0.5
  }, {
    name: "width",
    type: "number",
    default: 1
  }, {
    name: "height",
    type: "number",
    default: 1
  }],
  setImage: [{
    name: "url",
    type: "string",
    default: "'https://upload.wikimedia.org/wikipedia/en/thumb/8/80/Wikipedia-logo-v2.svg/1200px-Wikipedia-logo-v2.svg.png'"
  }, {
    name: "x",
    type: "number",
    default: -0.5
  }, {
    name: "y",
    type: "number",
    default: -0.5
  }, {
    name: "width",
    type: "number",
    default: 1
  }, {
    name: "height",
    type: "number",
    default: 1
  }],
  //2D transforms
  reset: [],
  arcsinh: [],
  arctanh: [],
  bent: [{
      name: "real",
      type: "number",
      default: 1
    },
    {
      name: "imaginary",
      type: "number",
      default: 1
    }
  ],
  blurCircle: [],
  blurGasket: [],
  blurGaussian: [{
    name: "pow",
    type: "number",
    default: 1
  }],
  blurNgon: [{
    name: "n",
    type: "number",
    default: "3"
  }],
  blurSine: [{
    name: "pow",
    type: "number",
    default: 1
  }],
  blurSquare: [],
  blurTriangle: [{
      name: "a",
      type: "complex",
      default: "C(-0.5, -Math.sqrt(3)*0.5)",
    },
    {
      name: "b",
      type: "complex",
      default: "C(-0.5, Math.sqrt(3)*0.5)",
    },
    {
      name: "c",
      type: "complex",
      default: "C(1, 0)",
    }
  ],
  bTransform: [{
      name: "rotate",
      type: "number",
      default: 0
    },
    {
      name: "power",
      type: "number",
      default: 1
    },
    {
      name: "move",
      type: "number",
      default: 0
    },
    {
      name: "split",
      type: "number",
      default: 0
    }
  ],
  bubble: [],
  circleInv: [],
  cpow: [{
      name: "rotation",
      type: "number",
      default: 1
    },
    {
      name: "depth",
      type: "number",
      default: 1
    },
    {
      name: "divisor",
      type: "number",
      default: 1
    },
    {
      name: "spread",
      type: "number",
      default: 1
    }
  ],
  cylinder: [],
  dc_poincareDisc: [{
      name: "p",
      type: "number",
      default: 4
    },
    {
      name: "q",
      type: "number",
      default: 5
    },
    {
      name: "iterations",
      type: "number",
      default: 10
    },
    {
      name: "checks0",
      type: "number",
      default: 1
    },
    {
      name: "checks1",
      type: "number",
      default: 1
    },
    {
      name: "checks2",
      type: "number",
      default: 1
    }
  ],
  disc: [],
  dragon: [{
      name: "a",
      type: "number",
      default: 3
    },
    {
      name: "b",
      type: "number",
      default: 3
    },
    {
      name: "c",
      type: "number",
      default: 3
    },
    {
      name: "bc",
      type: "number",
      default: 3
    },
    {
      name: "multiplier",
      type: "number",
      default: 1
    },
    {
      name: "horizontal",
      type: "number",
      default: 0.25
    },
    {
      name: "vertical",
      type: "number",
      default: 0.25
    },
    {
      name: "radial",
      type: "number",
      default: 0.5
    }
  ],
  ePush: [{
      name: "push",
      type: "number",
      default: 0
    },
    {
      name: "rotation",
      type: "number",
      default: 0
    }
  ],
  eRotate: [{
    name: "rotationAngle",
    type: "number",
    default: 0
  }],
  flipX: [],
  flipY: [],
  hypershape: [{
    name: "pow",
    type: "number",
    default: 3
  }],
  hypershift: [{
    name: "p",
    type: "complex",
    default: "C(0, 0)"
  }],
  hypertile3: [{
      name: "p",
      type: "number",
      default: 3
    },
    {
      name: "q",
      type: "number",
      default: 3
    },
    {
      name: "r",
      type: "number",
      default: 4
    },
    {
      name: "d",
      type: "number",
      default: 0.5
    }
  ],
  jac_cd: [{
    name: "k",
    type: "number",
    default: 0.5
  }],
  jac_cn: [{
    name: "k",
    type: "number",
    default: 0.5
  }],
  jac_dn: [{
    name: "k",
    type: "number",
    default: 0.5
  }],
  jac_elk: [{
    name: "k",
    type: "number",
    default: 0.5
  }],
  jac_sn: [{
    name: "k",
    type: "number",
    default: 0.5
  }],
  julian: [{
      name: "pow",
      type: "number",
      default: 1
    },
    {
      name: "dist",
      type: "number",
      default: 1
    }
  ],
  juliaq: [{
      name: "pow",
      type: "number",
      default: 1
    },
    {
      name: "div",
      type: "number",
      default: 1
    }
  ],
  juliascope: [{
      name: "pow",
      type: "number",
      default: 1
    },
    {
      name: "dist",
      type: "number",
      default: 1
    }
  ],
  mobius: [{
      name: "a",
      type: "complex",
      default: "C(1, 0)"
    },
    {
      name: "b",
      type: "complex",
      default: "C(0, 0)"
    },
    {
      name: "c",
      type: "complex",
      default: "C(0, 0)"
    },
    {
      name: "d",
      type: "complex",
      default: "C(1, 0)"
    }
  ],
  multiMobius: [{
      name: "a",
      type: "complex",
      default: "C(1, 0)"
    },
    {
      name: "b",
      type: "complex",
      default: "C(0, 0)"
    },
    {
      name: "c",
      type: "complex",
      default: "C(0, 0)"
    },
    {
      name: "d",
      type: "complex",
      default: "C(1, 0)"
    },
    {
      name: "iterations",
      type: "number",
      default: 1
    }
  ],
  murl2: [{
      name: "c",
      type: "number",
      default: 0
    },
    {
      name: "pow",
      type: "number",
      default: 2
    }
  ],
  nSplit: [{
      name: "n",
      type: "number",
      default: 3
    },
    {
      name: "split",
      type: "number",
      default: 0
    },
    {
      name: "wedge",
      type: "number",
      default: 1
    }
  ],
  pdj: [{
      name: "a",
      type: "number",
      default: 1
    },
    {
      name: "b",
      type: "number",
      default: 1
    },
    {
      name: "c",
      type: "number",
      default: 1
    },
    {
      name: "d",
      type: "number",
      default: 1
    }
  ],
  pointSymmetry: [{
      name: "centerX",
      type: "number",
      default: 0
    },
    {
      name: "centerY",
      type: "number",
      default: 0
    },
    {
      name: "order",
      type: "number",
      default: 1
    }
  ],
  rotate: [{
    name: "theta",
    type: "number",
    default: 0
  }],
  scale: [{
    name: "factor",
    type: "number",
    default: 1
  }],
  scale2: [{
    name: "real",
    type: "number",
    default: 1
  }, {
    name: "imaginary",
    type: "number",
    default: 1
  }],
  schwarzChristoffelMap: [{
      name: "n",
      type: "number",
      default: 3
    },
    {
      name: "k",
      type: "number",
      default: 0.5
    }
  ],
  schwarzChristoffelInverseMap: [{
      name: "n",
      type: "number",
      default: 3
    }
  ],
  schwarzTriangle: [{
      name: "alpha",
      type: "number",
      default: 60
    },
    {
      name: "beta",
      type: "number",
      default: 60
    },
    {
      name: "gamma",
      type: "number",
      default: 60
    }
  ],
  sinusoidal: [],
  skew: [{
      name: "real",
      type: "number",
      default: 0
    },
    {
      name: "imaginary",
      type: "number",
      default: 0
    }
  ],
  smartcrop: [{
      name: "power",
      type: "number",
      default: 4
    },
    {
      name: "radius",
      type: "number",
      default: 1
    },
    {
      name: "roundstr",
      type: "number",
      default: 0
    },
    {
      name: "roundwidth",
      type: "number",
      default: 1
    },
    {
      name: "distortion",
      type: "number",
      default: 1
    },
    {
      name: "cropmode",
      type: "number",
      default: 2
    }
  ],
  smartshape: [{
      name: "power",
      type: "number",
      default: 4
    },
    {
      name: "roundstr",
      type: "number",
      default: 0
    },
    {
      name: "roundwidth",
      type: "number",
      default: 1
    },
    {
      name: "distortion",
      type: "number",
      default: 1
    },
    {
      name: "compensation",
      type: "number",
      default: 1
    }
  ],
  splits: [{
      name: "real",
      type: "number",
      default: 0
    },
    {
      name: "imaginary",
      type: "number",
      default: 0
    }
  ],
  tileHelp: [],
  tileLog: [{
    name: "spread",
    type: "number",
    default: 1
  }],
  translate: [{
      name: "real",
      type: "number",
      default: 0
    },
    {
      name: "imaginary",
      type: "number",
      default: 0
    }
  ],
  trigCosh: [],
  trigExp: [],
  trigLog: [],
  trigSinh: [],
  trigTanh: [],
  unbubble: [],
  weierstrassElliptic: [{
      name: "w2",
      type: "complex",
      default: "C(0, 1)"
  }],
  //color transforms
  brighten: [{
    name: "amount",
    type: "number",
    default: 1
  }],
  color: [{
    name: "color",
    type: "object",
    default: {
      h: 0,
      s: 1,
      l: 0.5
    }
  }],
  gradient: [{
    name: "colorA",
    type: "array",
    default: [{
      red: 1,
      green: 1,
      blue: 1
    }, {
      red: 1,
      green: 0,
      blue: 0
    }]
  }],
  repeatingGradient: [{
    name: "colorA",
    type: "array",
    default: [{
      red: 1,
      green: 1,
      blue: 1
    }, {
      red: 1,
      green: 0,
      blue: 0
    }]
  }],
  hslShift: [{
    name: "h",
    type: "number",
    default: 0
  }, {
    name: "s",
    type: "number",
    default: 0
  }, {
    name: "l",
    type: "number",
    default: 0
  }],
  lerpHSL: [{
      name: "h",
      type: "number",
      default: 0
    },
    {
      name: "s",
      type: "number",
      default: 1
    },
    {
      name: "l",
      type: "number",
      default: 0.5
    }, {
      name: "weight",
      type: "number",
      default: 0.5
    }
  ],
  lerpColor: [{
    name: "color",
    type: "object",
    default: {
      red: 1,
      green: 1,
      blue: 1
    }
  }, {
    name: "weight",
    type: "number",
    default: 0.5
  }],
  lerpRGB: [{
      name: "r",
      type: "number",
      default: 1
    },
    {
      name: "g",
      type: "number",
      default: 1
    },
    {
      name: "b",
      type: "number",
      default: 1
    }, {
      name: "weight",
      type: "number",
      default: 0.5
    }
  ],
  normalizeColor: [],
  setHue: [{
    name: "hue",
    type: "number",
    default: 0
  }],
  setAlpha: [{
    name: "alpha",
    type: "number",
    default: 255
  }],
  setSaturation: [{
    name: "saturation",
    type: "number",
    default: 1
  }],
  setLightness: [{
    name: "lightness",
    type: "number",
    default: 0.5
  }],
};

const BUILT_IN_TRANSFORMS_DESCRIPTIONS = {
  identity: 'The identity function returns its input.',
  draw: 'Draws immediately. Useful for drawing multiple points per loop',
  dither: '(Shader only) Applies a Bayer dither to the image.',
  //3D models
  stl: 'renders a 3D model from an STL file encoded as a Base64 string',
  blurTetrahedron: 'Makes a unit radius tetrahedron',
  blurNormalCube: 'Makes a unit radius cube',
  blurOctahedron: 'Makes a unit radius octahedron',
  blurIcosahedron: 'Makes a unit radius icosahedron',
  blurDodecahedron: 'Makes a unit radius dodecahedron',
  blurTeapot: 'Makes a teapot',
  torus: 'Makes a 3D torus',
  umbilicTorus: 'TODO',
  //shaders
  shaderPass: '(Shader only) Calls the function passed for all pixels.',
  normalizeColors: '(Shader only) Normalizes the colors based on the brightest pixel',
  rainbowCirc: 'TODO',
  rainbowCircAdd: 'TODO',
  paletteMod: 'TODO',
  gamma: 'Applies a gamma to the color',
  ambientOcclusion: '(Shader only) gets the ambientOcclusion of a 3D scene',
  ambientOcclusion2: '(Shader only) gets the ambientOcclusion of a 3D scene',
  ambientOcclusionBig: '(Shader only) gets the ambientOcclusion of a 3D scene',
  edgeDetection: '(Shader only) shades the images based on edge detection',
  specular: '(Shader only) gets the specular lighting of a 3D scene',
  specularOrth: '(Shader only) gets the specular lighting of an orthographic 3D scene',
  normalMap: '(Shader only) gets a normal map of a 3D scene',
  matcap: '(Shader only) shades the scene using a given matcap image based on surface normals',
  heightMap: '(Shader only) gets a height map of a 3D scene',
  basicLighting: '(Shader only) applies very basic lighting to a 3D scene',
  basicEnvironment: '(Shader only) shades a 3D scene with basic enviormental lighting',
  basicEnvironmentOrth: '(Shader only) shades an orthographic 3D scene with basic enviormental lighting',
  advancedLighting: '(Shader only) shades a 3D scene with basic enviormental lighting',
  advancedLightingOrth: '(Shader only) shades an orthographic 3D scene with basic enviormental lighting',
  mist: '(Shader only) applies a volumetric mist effect to a 3D scene',
  //3D transforms
  mobius3D: 'TODO',
  hypershift3D: 'TODO',
  bubble3D: 'TODO',
  julian3D: 'TODO',
  juliaq3D: 'TODO',
  parabola3D: 'TODO',
  sphereInv: 'TODO',
  trigCosh3D: 'TODO',
  trigExp3D: 'TODO',
  trigLog3D: 'TODO',
  trigSinh3D: 'TODO',
  trigTanh3D: 'TODO',
  unbubble3D: 'TODO',
  matrix3D: 'TODO',
  blurSphere: 'Makes a unit radius sphere',
  blurSphereVolume: 'Makes a unit radius filled sphere',
  blurCube: 'Makes a unit width cube',
  blurCubeVolume: 'Makes a unit width filled cube',
  scale3D: 'scales uniformly in 3D',
  scale3D3: 'scales in 3D',
  translate3D: 'translates in 3D',
  rotate3D: 'rotates in 3D',
  perspective3D: 'applies a perspective transformation',
  viewBox: 'TODO',
  viewSphere: 'TODO',
  //Images
  blurImage: 'draws an image from a url',
  setImage: 'sets color based on an image from a url',
  //2D transforms
  reset: 'resets the pointer',
  arcsinh: 'TODO',
  arctanh: 'TODO',
  bent: 'TODO',
  blurCircle: 'draws a unit radius circle',
  blurGasket: 'TODO',
  blurGaussian: 'TODO',
  blurNgon: 'draws a unit radius N-gon',
  blurSine: 'TODO',
  blurSquare: 'draws a unit width square',
  blurTriangle: 'draws a triangle',
  bTransform: 'TODO',
  bubble: 'TODO',
  circleInv: 'TODO',
  cpow: 'TODO',
  cylinder: 'TODO',
  dc_poincareDisc: 'TODO',
  disc: 'TODO',
  dragon: 'TODO',
  ePush: 'TODO',
  eRotate: 'TODO',
  flipX: 'TODO',
  flipY: 'TODO',
  hypershape: 'TODO',
  hypershift: 'TODO',
  hypertile3: 'TODO',
  jac_cd: 'TODO',
  jac_cn: 'TODO',
  jac_dn: 'TODO',
  jac_elk: 'TODO',
  jac_sn: 'TODO',
  julian: 'TODO',
  juliaq: 'TODO',
  juliascope: 'TODO',
  mobius: 'TODO',
  multiMobius: 'TODO',
  murl2: 'TODO',
  nSplit: 'TODO',
  pdj: 'TODO',
  pointSymmetry: 'TODO',
  rotate: 'rotates in 2D',
  scale: 'scales uniformly in 2D',
  scale2: 'scales in 2D',
  schwarzChristoffelMap: 'TODO',
  schwarzChristoffelInverseMap: 'TODO',
  schwarzTriangle: 'TODO',
  sinusoidal: 'TODO',
  skew: 'TODO',
  smartcrop: 'TODO',
  smartshape: 'TODO',
  splits: 'TODO',
  tileHelp: 'TODO',
  tileLog: 'TODO',
  translate: 'translates in 2D',
  trigCosh: 'TODO',
  trigExp: 'TODO',
  trigLog: 'TODO',
  trigSinh: 'TODO',
  trigTanh: 'TODO',
  unbubble: 'TODO',
  weierstrassElliptic: 'TODO',
  //color transfomrs
  brighten: 'TODO',
  color: 'sets the color',
  gradient: 'TODO',
  repeatingGradient: 'TODO',
  hslShift: 'shifts the Hue/Saturation/Lightness by a given amount',
  lerpColor: 'blends towards a given color',
  lerpHSL: 'blends towards a given color in HSL',
  lerpRGB: 'blends towards a given color in RGB',
  normalizeColor: 'normalizes color to be maximally bright',
  setHue: 'sets the hue',
  setAlpha: 'sets the alpha',
  setSaturation: 'sets the saturation',
  setLightness: 'sets the lightness',
};

const BUILT_IN_TRANSFORMS_PREVIEW_CODE = {
  identity: {
    type: 'any',
    default: 'Checkered'
  },
  draw: {
    type: 'normal',
    default: 'Checkered'
  },
  dither: {
    type: 'shader',
    default: 'Dark Gray Shader'
  },
  //3D models
  stl: {
    type: 'normal',
    default: '3D Teapot'
  },
  blurTetrahedron: {
    type: 'normal',
    default: '3D Teapot'
  },
  blurNormalCube: {
    type: 'normal',
    default: '3D Teapot'
  },
  blurOctahedron: {
    type: 'normal',
    default: '3D Teapot'
  },
  blurIcosahedron: {
    type: 'normal',
    default: '3D Teapot'
  },
  blurDodecahedron: {
    type: 'normal',
    default: '3D Teapot'
  },
  blurTeapot: {
    type: 'normal',
    default: '3D Teapot'
  },
  torus: {
    type: 'normal',
    default: '3D Teapot'
  },
  umbilicTorus: {
    type: 'normal',
    default: '3D Teapot'
  },
  //shaders
  shaderPass: {
    type: 'shader',
    default: 'Checkered Shader'
  },
  normalizeColors: {
    type: 'shader',
    default: 'Dark Gray Shader'
  },
  rainbowCirc: {
    type: 'normal',
    default: 'Checkered'
  },
  rainbowCircAdd: {
    type: 'normal',
    default: 'Checkered'
  },
  paletteMod: {
    type: 'normal',
    default: 'Checkered'
  },
  gamma: {
    type: 'shader',
    default: 'Checkered Shader'
  },
  ambientOcclusion: {
    type: 'shader',
    default: '3D Teapot Shader'
  },
  ambientOcclusion2: {
    type: 'shader',
    default: '3D Teapot Shader'
  },
  ambientOcclusionBig: {
    type: 'shader',
    default: '3D Teapot Shader'
  },
  edgeDetection: {
    type: 'shader',
    default: '3D Teapot Shader'
  },
  specular: {
    type: 'shader',
    default: '3D Teapot Shader'
  },
  specularOrth: {
    type: 'shader',
    default: '3D Teapot Orth Shader'
  },
  normalMap: {
    type: 'shader',
    default: '3D Teapot Orth Shader'
  },
  matcap: {
    type: 'shader',
    default: '3D Teapot Shader'
  },
  heightMap: {
    type: 'shader',
    default: '3D Teapot Orth Shader'
  },
  basicLighting: {
    type: 'shader',
    default: '3D Teapot Shader'
  },
  basicEnvironment: {
    type: 'shader',
    default: '3D Teapot Shader'
  },
  basicEnvironmentOrth: {
    type: 'shader',
    default: '3D Teapot Orth Shader'
  },
  advancedLighting: {
    type: 'shader',
    default: '3D Teapot Shader'
  },
  advancedLightingOrth: {
    type: 'shader',
    default: '3D Teapot Orth Shader'
  },
  mist: {
    type: 'shader',
    default: '3D Teapot Shader'
  },
  //3D transforms
  mobius3D: {
    type: 'normal',
    default: '3D Teapot'
  },
  hypershift3D: {
    type: 'normal',
    default: '3D Teapot'
  },
  bubble3D: {
    type: 'normal',
    default: '3D Teapot'
  },
  julian3D: {
    type: 'normal',
    default: '3D Teapot'
  },
  juliaq3D: {
    type: 'normal',
    default: '3D Teapot'
  },
  parabola3D: {
    type: 'normal',
    default: '3D Teapot'
  },
  sphereInv: {
    type: 'normal',
    default: '3D Teapot'
  },
  trigCosh3D: {
    type: 'normal',
    default: '3D Teapot'
  },
  trigExp3D: {
    type: 'normal',
    default: '3D Teapot'
  },
  trigLog3D: {
    type: 'normal',
    default: '3D Teapot'
  },
  trigSinh3D: {
    type: 'normal',
    default: '3D Teapot'
  },
  trigTanh3D: {
    type: 'normal',
    default: '3D Teapot'
  },
  unbubble3D: {
    type: 'normal',
    default: '3D Teapot'
  },
  matrix3D: {
    type: 'normal',
    default: '3D Teapot'
  },
  blurSphere: {
    type: 'normal',
    default: '3D Teapot'
  },
  blurSphereVolume: {
    type: 'normal',
    default: '3D Teapot'
  },
  blurCube: {
    type: 'normal',
    default: '3D Teapot'
  },
  blurCubeVolume: {
    type: 'normal',
    default: '3D Teapot'
  },
  scale3D: {
    type: 'normal',
    default: '3D Teapot'
  },
  scale3D3: {
    type: 'normal',
    default: '3D Teapot'
  },
  translate3D: {
    type: 'normal',
    default: '3D Teapot'
  },
  rotate3D: {
    type: 'normal',
    default: '3D Teapot'
  },
  perspective3D: {
    type: 'normal',
    default: '3D Teapot'
  },
  viewBox: {
    type: 'normal',
    default: '3D Teapot'
  },
  viewSphere: {
    type: 'normal',
    default: '3D Teapot'
  },
  //Images
  blurImage: {
    type: 'normal',
    default: 'Checkered'
  },
  setImage: {
    type: 'normal',
    default: 'Checkered'
  },
  //2D transforms
  reset: {
    type: 'normal',
    default: 'Checkered'
  },
  arcsinh: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  arctanh: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  bent: {
    type: 'normal',
    default: 'Checkered'
  },
  blurCircle: {
    type: 'normal',
    default: 'Checkered'
  },
  blurGasket: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  blurGaussian: {
    type: 'normal',
    default: 'Checkered'
  },
  blurNgon: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  blurSine: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  blurSquare: {
    type: 'normal',
    default: 'Checkered'
  },
  blurTriangle: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  bTransform: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  bubble: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  circleInv: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  cpow: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  cylinder: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  dc_poincareDisc: {
    type: 'normal',
    default: 'Checkered'
  },
  disc: {
    type: 'normal',
    default: 'Checkered'
  },
  dragon: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  ePush: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  eRotate: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  flipX: {
    type: 'normal',
    default: '3D Teapot'
  },
  flipY: {
    type: 'normal',
    default: '3D Teapot'
  },
  hypershape: {
    type: 'normal',
    default: 'Checkered'
  },
  hypershift: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  hypertile3: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  jac_cd: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  jac_cn: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  jac_dn: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  jac_elk: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  jac_sn: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  julian: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  juliaq: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  juliascope: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  mobius: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  multiMobius: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  murl2: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  nSplit: {
    type: 'normal',
    default: 'Checkered'
  },
  pdj: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  pointSymmetry: {
    type: 'normal',
    default: 'Checkered'
  },
  rotate: {
    type: 'normal',
    default: 'Checkered'
  },
  scale: {
    type: 'normal',
    default: 'Checkered'
  },
  scale2: {
    type: 'normal',
    default: 'Checkered'
  },
  schwarzChristoffelMap: {
    type: 'normal',
    default: 'Checkered'
  },
  schwarzChristoffelInverseMap: {
    type: 'normal',
    default: 'Checkered'
  },
  schwarzTriangle: {
    type: 'normal',
    default: 'Checkered'
  },
  sinusoidal: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  skew: {
    type: 'normal',
    default: 'Checkered'
  },
  smartcrop: {
    type: 'normal',
    default: 'Checkered'
  },
  smartshape: {
    type: 'normal',
    default: 'Checkered'
  },
  splits: {
    type: 'normal',
    default: 'Checkered'
  },
  tileHelp: {
    type: 'normal',
    default: 'Checkered'
  },
  tileLog: {
    type: 'normal',
    default: 'Checkered'
  },
  translate: {
    type: 'normal',
    default: 'Checkered'
  },
  trigCosh: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  trigExp: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  trigLog: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  trigSinh: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  trigTanh: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  unbubble: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  weierstrassElliptic: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  //color transfomrs
  brighten: {
    type: 'normal',
    default: 'Checkered'
  },
  color: {
    type: 'normal',
    default: 'Checkered'
  },
  gradient: {
    type: 'normal',
    default: 'Checkered'
  },
  repeatingGradient: {
    type: 'normal',
    default: 'Checkered Plane'
  },
  hslShift: {
    type: 'normal',
    default: 'Checkered'
  },
  lerpColor: {
    type: 'normal',
    default: 'Checkered'
  },
  lerpHSL: {
    type: 'normal',
    default: 'Checkered'
  },
  lerpRGB: {
    type: 'normal',
    default: 'Checkered'
  },
  normalizeColor: {
    type: 'normal',
    default: 'Checkered'
  },
  setHue: {
    type: 'normal',
    default: 'Checkered'
  },
  setAlpha: {
    type: 'normal',
    default: 'Checkered'
  },
  setSaturation: {
    type: 'normal',
    default: 'Checkered'
  },
  setLightness: {
    type: 'normal',
    default: 'Checkered'
  },
};
