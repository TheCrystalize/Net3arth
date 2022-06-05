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
  let total = oldBuffer.z + 1;
  if(newBuffer.z > 0){
    add = oldBuffer.z + newBuffer.z;
  }
  return {
    red: (oldBuffer.red * oldBuffer.z + newBuffer.red) / total,
    green: (oldBuffer.green * oldBuffer.z + newBuffer.green) / total,
    blue: (oldBuffer.blue * oldBuffer.z + newBuffer.blue) / total,
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
    if(z.re === z.width-1 && z.im === z.height-1) {
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
      red: newBuffer[0][(z.re - 1) + z.im * z.width],
      green: newBuffer[1][(z.re - 1) + z.im * z.width],
      blue: newBuffer[2][(z.re - 1) + z.im * z.width],
      z: newZBuffer[(z.re - 1) + z.im * z.width],
    };
  }
}

function normalizeColors() {
  let brightness;
  return z => {
    if(z.re === z.width-1 && z.im === z.height-1) {
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
      green:  z.green / brightness,
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

function pow(z, p) {
  let n = p * Math.atan2(z.im, z.re);
  let m = Math.exp(p * Math.log(z.re * z.re + z.im * z.im));
  return {
    re: Math.cos(n) * m,
    im: Math.sin(n) * m
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

function tcPow(z, c) {
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
    prod = mult(prod, tcPow(C(1 + 1 / i, 0), z));
    prod = div(prod, addScalar(div(z, C(i, 0)), 1))
  }
  gammaMap.set(z, prod)
  return prod
}

function cGamma(z) {
  let eps = 0.01;
  let sum = C(0, 0);
  for(let i = 0; i < 10; i += eps) {
    sum = add(sum, multScalar(mult(tcPow(C(i, 0), addScalar(z, -1)), exp(C(-i, 0))), eps));
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

function round(x) {
  let f = Math.floor(x)
  let n = x > 0 ? (Math.abs(x - f) > 0.5 ? 1 : 0) : (Math.abs(x - f) > 0.5 ? -1 : 0);
  return f + n;
}

function ccpow(z, p) {
  return exp(multScalar(log(z), p))
}

function hypergeo2F1PowSerI(a, b, c, z) {
  let y = C(0, 0);
  let yp = C(0, 0);
  for(let k = 0; k < 40; k++) {
    yp = divScalar(multScalar(ccpow(z, k), pochhammer(a, k) * pochhammer(b, k)), gammaLanczos(c + k) * factorial(k));
    y = add(y, yp);
  }
  return y;
}

function powSerO(a, b, c, z) {
  let y = C(0, 0);
  let yp = C(0, 0);
  for(let k = 0; k < 40; k++) {
    yp = divScalar(multScalar(ccpow(z, -k), pochhammer(a, k) * pochhammer(a - c + 1, k)), factorial(k) * gammaLanczos(a - b + k + 1));
    y = add(y, yp);
  }
  return y;
}

function hypergeo2F1PowSerO(a, b, c, z) {
  let sc = Math.PI / Math.sin(Math.PI * (b - a));
  let px = divScalar(ccpow(neg(z), -a), gammaLanczos(b) * gammaLanczos(c - a));
  let py = divScalar(ccpow(neg(z), -b), gammaLanczos(a) * gammaLanczos(c - b));
  let x = powSerO(a, b, c, z);
  let y = powSerO(b, a, c, z);
  return multScalar(add(mult(px, x), neg(mult(py, y))), sc)
}

function hypergeo2F1(a, b, c, z) {
  let y = C(0, 0);
  if(Math.hypot(z.re, z.im) > 1) {
    y = hypergeo2F1PowSerO(a, b, c, z);
  } else {
    y = hypergeo2F1PowSerI(a, b, c, z);
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

function scHGt(z, n, k) {
  let nums = [0, 0,
    2 + n, (1 + n) * (2 + n), (1 + n) * (2 + n) * (2 + 3 * n), (1 + n) * (2 + n) * (1 + 2 * n) * (2 + 3 * n), (1 + n) * (2 + n) * (1 + 2 * n) * (2 + 3 * n) * (2 + 5 * n),
    (1 + n) * (2 + n) * (1 + 2 * n) * (1 + 3 * n) * (2 + 3 * n) * (2 + 5 * n), (1 + n) * (2 + n) * (1 + 2 * n) * (1 + 3 * n) * (2 + 3 * n) * (2 + 5 * n) * (2 + 7 * n), (1 + n) * (2 + n) * (1 + 2 * n) * (1 + 3 * n) * (1 + 4 * n) * (2 + 3 * n) * (2 + 5 * n) * (2 + 7 * n), (1 + n) * (2 + n) * (1 + 2 * n) * (1 + 3 * n) * (1 + 4 * n) * (2 + 3 * n) * (2 + 5 * n) * (2 + 7 * n) * (2 + 9 * n), (1 + n) * (2 + n) * (1 + 2 * n) * (1 + 3 * n) * (1 + 4 * n) * (2 + 3 * n) * (2 + 5 * n) * (2 + 7 * n) * (2 + 9 * n),
    (1 + n) * (2 + n) * (1 + 2 * n) * (1 + 3 * n) * (1 + 4 * n) * (2 + 3 * n) * (2 + 5 * n) * (2 + 7 * n) * (2 + 9 * n) * (2 + 11 * n), (1 + n) * (2 + n) * (1 + 2 * n) * (1 + 3 * n) * (1 + 4 * n) * (2 + 3 * n) * (2 + 5 * n) * (2 + 7 * n) * (2 + 9 * n) * (2 + 11 * n), (1 + n) * (2 + n) * (1 + 2 * n) * (1 + 3 * n) * (1 + 4 * n) * (2 + 3 * n) * (2 + 5 * n) * (2 + 7 * n) * (2 + 9 * n) * (2 + 11 * n) * (2 + 13 * n), (1 + n) * (2 + n) * (1 + 2 * n) * (1 + 3 * n) * (1 + 4 * n) * (2 + 3 * n) * (2 + 5 * n) * (2 + 7 * n) * (2 + 9 * n) * (2 + 11 * n) * (2 + 13 * n), (1 + n) * (2 + n) * (1 + 2 * n) * (1 + 3 * n) * (1 + 4 * n) * (2 + 3 * n) * (2 + 5 * n) * (2 + 7 * n) * (2 + 9 * n) * (2 + 11 * n) * (2 + 13 * n) * (2 + 15 * n),
    (1 + n) * (2 + n) * (1 + 2 * n) * (1 + 3 * n) * (1 + 4 * n) * (2 + 3 * n) * (2 + 5 * n) * (2 + 7 * n) * (2 + 9 * n) * (2 + 11 * n) * (2 + 13 * n) * (2 + 15 * n), (1 + n) * (2 + n) * (1 + 2 * n) * (1 + 3 * n) * (1 + 4 * n) * (2 + 3 * n) * (2 + 5 * n) * (2 + 7 * n) * (2 + 9 * n) * (2 + 11 * n) * (2 + 13 * n) * (2 + 15 * n) * (2 + 17 * n), (1 + n) * (2 + n) * (1 + 2 * n) * (1 + 3 * n) * (1 + 4 * n) * (2 + 3 * n) * (2 + 5 * n) * (2 + 7 * n) * (2 + 9 * n) * (2 + 11 * n) * (2 + 13 * n) * (2 + 15 * n) * (2 + 17 * n), (1 + n) * (2 + n) * (1 + 2 * n) * (1 + 3 * n) * (1 + 4 * n) * (2 + 3 * n) * (2 + 5 * n) * (2 + 7 * n) * (2 + 9 * n) * (2 + 11 * n) * (2 + 13 * n) * (2 + 15 * n) * (2 + 17 * n) * (2 + 19 * n), (1 + n) * (2 + n) * (1 + 2 * n) * (1 + 3 * n) * (1 + 4 * n) * (2 + 3 * n) * (2 + 5 * n) * (2 + 7 * n) * (2 + 9 * n) * (2 + 11 * n) * (2 + 13 * n) * (2 + 15 * n) * (2 + 17 * n) * (2 + 19 * n)
  ];
  let dens = [0, 0,
    Math.pow(n, 2) * (1 + 2 * n), 3 * Math.pow(n, 3) * (1 + 3 * n), 6 * Math.pow(n, 4) * (1 + 4 * n), 15 * Math.pow(n, 5) * (1 + 5 * n), 90 * Math.pow(n, 6) * (1 + 6 * n),
    315 * Math.pow(n, 7) * (1 + 2 * n), 2520 * Math.pow(n, 8) * (1 + 8 * n), 11340 * Math.pow(n, 9) * (1 + 9 * n), 113400 * Math.pow(n, 10) * (1 + 10 * n), 623700 * Math.pow(n, 11) * (1 + 11 * n),
    7484400 * Math.pow(n, 12) * (1 + 12 * n), 48648600 * Math.pow(n, 13) * (1 + 13 * n), 681080400 * Math.pow(n, 14) * (1 + 14 * n), 5108103000 * Math.pow(n, 15) * (1 + 15 * n), 81729648000 * Math.pow(n, 16) * (1 + 16 * n),
    694702008000 * Math.pow(n, 17) * (1 + 17 * n), 12504636144000 * Math.pow(n, 18) * (1 + 18 * n), 118794043368000 * Math.pow(n, 19) * (1 + 19 * n), 2375880867360000 * Math.pow(n, 20) * (1 + 20 * n), 24946749107280000 * Math.pow(n, 21) * (1 + 21 * n)
  ];
  let sum = C(1, 0);
  let fin = div(multScalar(ccpow(z, n), 2), C(n + Math.pow(n, 2), 0));
  for(let i = 2; i < 22; i++) {
    sum = add(sum, div(multScalar(ccpow(z, i * n), (i == 3 ? 2 : 1) * nums[i]), C(dens[i], 0)));
  }
  return mult(z, add(sum, fin));
}

function hypergeo2F1Coefficients(a, b, c, nomial) {
  let poly = [nomial];
  for(let i = 0; i < nomial; i++) {
    poly[i] = C((dynamicPochhammer(a, i) * dynamicPochhammer(b, i)) / (dynamicPochhammer(c, i) * dynamicFactorial(i)), 0);
  }
  return poly
}

function hypergeo2F1pn(z, poly) {
  let sum = C(0, 0);
  for(let i = 0; i < poly.length; i++) {
    sum = add(sum, mult(ccpow(z, i), poly[i]))
  }
  return sum
}

function beta(x, y) {
  return div(mult(dynamicGamma(C(x, 0)), dynamicGamma(C(y, 0))), dynamicGamma(C(x + y, 0)))
}

function scFunc(a, b, c, x, z) {
  let fr = Math.pow(x, b - 1) * Math.pow(1 - x, c - b - 1);
  let cr = ccpow(addScalar(neg(multScalar(z, x)), 1), -a)
  return multScalar(cr, fr)
}

function integrateSCFunc(a, b, c, z) {
  let eps = 0.001;
  let sum = C(0, 0);
  for(let i = eps; i < 1; i += eps) {
    sum = add(sum, multScalar(scFunc(a, b, c, i, z), eps))
  }
  return sum
}

function hypergeoViaIntregration(a, b, c, z) {
  return div(integrateSCFunc(a, b, c, z), beta(b, c - b))
}

function hypergeo2F1CoefficientsCoord(a, b, c, z0, nomial) {
  let poly = [nomial];
  for(let i = 0; i < nomial; i++) {
    //let poly0 = hypergeo2F1Coefficients(a+i,b+i,c+i,nomial)
    poly[i] = C((dynamicPochhammer(a, i) * dynamicPochhammer(b, i)) / (dynamicPochhammer(c, i) * dynamicFactorial(i)), 0);
    //poly[i] *= hypergeo2F1pn(z0, poly0);
    poly[i] = mult(poly[i], hypergeoViaIntregration(a + i, b + i, c + i, z0))
  }
  return poly
}

function hypergeo2F1fast(z, polyMap) {
  let minDist = add(polyMap.keys().next().value, neg(z));
  let minDistance = Math.hypot(minDist.re, minDist.im);
  let minZ = polyMap.keys().next().value;
  for(const [key, value] of polyMap.entries()) {
    let keyMod = Math.hypot(key.re - z.re, key.im - z.im);
    if(keyMod < minDistance) {
      minDistance = keyMod;
      minZ = key;
    }
  }
  return hypergeo2F1pn(z, polyMap.get(minZ))
}

function scHG(z, n, k, poly) {
  return mult(z, hypergeo2F1pn(ccpow(z, n), poly))
}

function scHGfast(z, n, k, polyMap) {
  return mult(z, hypergeo2F1fast(ccpow(z, n), polyMap))
}

/* (c^(
    20 n) (1 + n) (2 + n) (1 + 2 n) (1 + 3 n) (2 + 3 n) (1 +
      4 n) (1 + 5 n) (2 + 5 n) (1 + 6 n) (1 + 7 n) (2 + 7 n) (1 +
      8 n) (1 + 9 n) (2 + 9 n) (2 + 11 n) (2 + 13 n) (2 + 15 n) (2 +
      17 n) (2 + 19 n))/(2375880867360000 n^20 (1 + 20 n))
*/

/*function SCinteg(x, y, func) {
  let eps = 0.1;

  let gammaPrime = C(x,y);
  let sum = C(0,0);
  for(let s = 0; s < 1; s += eps) {
    let first = mult(func(C(x * s, y * s)), multScalar(gammaPrime, eps));
    let second = mult(func(C(x * s+eps, y * s+eps)), multScalar(gammaPrime, eps))
    let third = mult(func(C(x * s+eps*2, y * s+eps*2)), multScalar(gammaPrime, eps));
    let fourth = mult(func(C(x * s+eps*3, y * s+eps*3)), multScalar(gammaPrime, eps));
    let fifth = mult(func(C(x * s+eps*4, y * s+eps*4)), multScalar(gammaPrime, eps));
    sum = add(sum, multScalar(add(add(add(add(multScalar(first,7),multScalar(second,32)),multScalar(third,12)),multScalar(fourth,32)),multScalar(fifth,7)),2/180));
  }
  return sum
}
/*transforms*/


function schwarzChristoffelmap(n, k) {
  let nom = 85;
  let dec = 0.5;
  let poly = hypergeo2F1Coefficients(1 / n, 2 / n, 1 + 1 / n, nom)
  let poly2 = hypergeo2F1CoefficientsCoord(1 / n, 2 / n, 1 + 1 / n, C(dec, 0), nom);
  let poly3 = hypergeo2F1CoefficientsCoord(1 / n, 2 / n, 1 + 1 / n, C(Math.sqrt(dec), Math.sqrt(dec)), nom);
  let poly4 = hypergeo2F1CoefficientsCoord(1 / n, 2 / n, 1 + 1 / n, C(0, dec), nom);
  let poly5 = hypergeo2F1CoefficientsCoord(1 / n, 2 / n, 1 + 1 / n, C(-Math.sqrt(dec), Math.sqrt(dec)), nom);
  let poly6 = hypergeo2F1CoefficientsCoord(1 / n, 2 / n, 1 + 1 / n, C(-dec, 0), nom);
  let poly7 = hypergeo2F1CoefficientsCoord(1 / n, 2 / n, 1 + 1 / n, C(-Math.sqrt(dec), -Math.sqrt(dec)), nom);
  let poly8 = hypergeo2F1CoefficientsCoord(1 / n, 2 / n, 1 + 1 / n, C(0, -dec), nom);
  let poly9 = hypergeo2F1CoefficientsCoord(1 / n, 2 / n, 1 + 1 / n, C(Math.sqrt(dec), -Math.sqrt(dec)), nom);
  let poly10 = hypergeo2F1CoefficientsCoord(1 / n, 2 / n, 1 + 1 / n, C(-0.5, Math.sqrt(3) * 0.5), nom);
  let poly11 = hypergeo2F1CoefficientsCoord(1 / n, 2 / n, 1 + 1 / n, C(-0.5, -Math.sqrt(3) * 0.5), nom);
  let poly12 = hypergeo2F1CoefficientsCoord(1 / n, 2 / n, 1 + 1 / n, C(Math.sqrt(3), 0), nom);
  let theMap = new Map();
  theMap.set(C(0, 0), poly);
  theMap.set(C(dec, 0), poly2);
  theMap.set(C(Math.sqrt(dec), Math.sqrt(dec)), poly3);
  theMap.set(C(0, dec), poly4);
  theMap.set(C(-Math.sqrt(dec), Math.sqrt(dec)), poly5);
  theMap.set(C(-dec, 0), poly6);
  theMap.set(C(-Math.sqrt(dec), -Math.sqrt(dec)), poly7);
  theMap.set(C(0, -dec), poly8);
  theMap.set(C(Math.sqrt(dec), -Math.sqrt(dec)), poly9);
  theMap.set(C(-0.5, Math.sqrt(3) * 0.5), poly10);
  theMap.set(C(-0.5, -Math.sqrt(3) * 0.5), poly11);
  theMap.set(C(Math.sqrt(3), 0), poly12);
  console.log(theMap)
  return z => {
    return {
      ...z,
      ...scHGfast(z, n, k, theMap)
    }
  }
}

/*function schwarzChristoffelmap(n, k) {
  let poly = hypergeo2F1Coefficients(1/n, 2/n, 1+1/n, 80)
  return z => {
    return {
      ...z,
      ...scHG(z,n,k,poly)
    }
  }
}


/*function schwarzChristoffelmap(n, k) {
  function funcyN(c) {
      return ccpow(add(C(1,0),neg(ccpow(c, n))),-2/n);
  }
  return z => {
    return {
      ...z,
      ...SCinteg(z.re, z.im, funcyN)
    }
  }
}
*/
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
  return z => {
    let a = Math.random() * Math.PI * 2,
      r = Math.sqrt(Math.random());
    return {
      ...z,
      re: Math.cos(a) * r,
      im: Math.sin(a) * r,
      z: 0,
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
  return z => {
    return {
      ...z,
      re: Math.random() - 0.5,
      im: Math.random() - 0.5,
      z: 0,
    }
  }
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

function jac_cn(k) {
  return z => {
    let jx = jacElliptic(z.re, k);
    let jy = jacElliptic(z.im, 1 - k);

    let numx = jx.c * jy.c;
    let numy = jx.d * jx.s * jy.d * jy.s;

    let denom = jx.s ** 2 * jy.s ** 2 * k + jy.c ** 2;
    denom = 1 / (denom);

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

    let denom = jx.s ** 2 * jy.s ** 2 * k + jy.c ** 2;
    denom = 1 / (denom);

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

    let denom = jx.s ** 2 * jy.s ** 2 * k + jy.c ** 2;
    denom = 1 / (denom);

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
    let nz = ccpow(z, alph);
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
  return n === 0 ? -Number.MAX_SAFE_INTEGER : n;
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
  return z => {
    let u = Math.random();
    let v = Math.random();
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
  return z => {
    let face = Math.random() * 6 >> 0;
    let x = Math.random() - 0.5;
    let y = Math.random() - 0.5;
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

function ambientOcclusion(s) {
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
        if(z.re + i < 0 || z.re + i >= z.width || z.im + j < 0 || z.im + j >= z.height){
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

function ambientOcclusionBig(s, step, depth) {
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
        if(z.re + i * step < 0 || z.re + i * step >= z.width || z.im + j * step < 0 || z.im + j * step >= z.height){
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

function reflect(vector, normal) {
  // v - 2 * (v dot n) * n
  return vectorSum(vector, vectorTimes(normal, -2 * (vector[0] * normal[0] + vector[1] * normal[1] + vector[2] * normal[2])));
}

function specularOrth(theta1, theta2, ior, environment, p) {
  let skyBoxLights = [lightRoomLights, dayLights][environment];
  let rotation1 = rotate3D(theta1, 0, 0);
  let rotation2 = rotate3D(0, theta2, 0);

  return z => {
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
  let skyBox = [lightRoomEnvironment, dayEnvironment][environment];
  let skyBoxLights = [lightRoomLights, dayLights][environment];
  let rotation1 = rotate3D(theta1, 0, 0);
  let rotation2 = rotate3D(0, theta2, 0);

  return z => {
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
  let skyBox = [lightRoomEnvironment, dayEnvironment][environment];
  let skyBoxLights = [lightRoomLights, dayLights][environment];
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
  let skyBoxLights = [lightRoomLights, dayLights][environment];
  let rotation1 = rotate3D(theta1, 0, 0);
  let rotation2 = rotate3D(0, theta2, 0);

  return z => {
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
  let skyBox = [lightRoomEnvironment, dayEnvironment][environment];
  let skyBoxLights = [lightRoomLights, dayLights][environment];
  let rotation1 = rotate3D(theta1, 0, 0);
  let rotation2 = rotate3D(0, theta2, 0);

  return z => {
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
  let skyBox = [lightRoomEnvironment, dayEnvironment][environment];
  let skyBoxLights = [lightRoomLights, dayLights][environment];
  let rotation1 = rotate3D(theta1, 0, 0);
  let rotation2 = rotate3D(0, theta2, 0);

  return z => {
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
  dither: dither,
  //shaders
  shaderPass: shaderPass,
  rainbowCirc: rainbowCirc,
  rainbowCircAdd: rainbowCircAdd,
  paletteMod: paletteMod,
  gamma: gamma,
  ambientOcclusion: ambientOcclusion,
  ambientOcclusionBig: ambientOcclusionBig,
  specular: specular,
  specularOrth: specularOrth,
  normalMap: normalMap,
  heightMap: heightMap,
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
  blurCube: blurCube,
  scale3D: scale3D,
  scale3D3: scale3D3,
  translate3D: translate3D,
  rotate3D: rotate3D,
  perspective3D: perspective3D,
  viewBox: viewBox,
  viewSphere: viewSphere,
  //2D transforms
  reset: reset,
  arcsinh: arcsinh,
  arctanh: arctanh,
  bent: bent,
  blurCircle: blurCircle,
  blurGasket: blurGasket,
  blurGaussian: blurGaussian,
  blurSine: blurSine,
  blurSquare: blurSquare,
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
  schwarzChristoffelmap: schwarzChristoffelmap,
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
  dither: [{
    name: "matrix size",
    type: "number",
    default: 5
  }],
  //shaders
  shaderPass: [{
    name: "transform",
    type: "function",
    default: "heightMap()"
  }],
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
      name: "_z",
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
      name: "_z",
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
    name: "sample size",
    type: "number",
    default: 5
  }],
  ambientOcclusionBig: [{
    name: "sample size",
    type: "number",
    default: 5
  },{
    name: "step size",
    type: "number",
    default: 2
  },{
    name: "depth",
    type: "number",
    default: 0
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
    default: 0,
  }, {
    name: "halfLength",
    type: "number",
    default: 0.5,
  }, {
    name: "mistColor",
    type: "object",
    default: {
      red: 0.7,
      green: 0.7,
      blue: 0.8
    },
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
      name: "_z",
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
  blurCube: [],
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
  blurSine: [{
    name: "pow",
    type: "number",
    default: 1
  }],
  blurSquare: [],
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
    default: {
      re: 0,
      im: 0
    }
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
      default: {
        re: 1,
        im: 0
      }
    },
    {
      name: "b",
      type: "complex",
      default: {
        re: 0,
        im: 0
      }
    },
    {
      name: "c",
      type: "complex",
      default: {
        re: 0,
        im: 0
      }
    },
    {
      name: "d",
      type: "complex",
      default: {
        re: 1,
        im: 0
      }
    }
  ],
  multiMobius: [{
      name: "a",
      type: "complex",
      default: {
        re: 1,
        im: 0
      }
    },
    {
      name: "b",
      type: "complex",
      default: {
        re: 0,
        im: 0
      }
    },
    {
      name: "c",
      type: "complex",
      default: {
        re: 0,
        im: 0
      }
    },
    {
      name: "d",
      type: "complex",
      default: {
        re: 1,
        im: 0
      }
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
  schwarzChristoffelmap: [{
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
  schwarzTriangle: [{
      name: "_alph",
      type: "number",
      default: 60
    },
    {
      name: "_bet",
      type: "number",
      default: 60
    },
    {
      name: "_gam",
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
