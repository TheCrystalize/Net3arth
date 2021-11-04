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

function zBuffer(oldBuffer, newBuffer) {
  if(oldBuffer.z !== 0 && oldBuffer.z > newBuffer.z) {
    return oldBuffer;
  }
  return {
    red: newBuffer.red * newBuffer.alpha,
    green: newBuffer.green * newBuffer.alpha,
    blue: newBuffer.blue * newBuffer.alpha,
    z: newBuffer.z
  };
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
  const n = p * Math.atan2(z.im, z.re);
  return {
    re: Math.cos(n) * Math.exp(p * Math.log(z.re * z.re + z.im * z.im)),
    im: Math.sin(n) * Math.exp(p * Math.log(z.re * z.re + z.im * z.im))
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

/*transforms*/

function reset() {
  return z => {
    return {
      re: 0.000001,
      im: 0.000002,
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
    let r = (pow == 1 ? Math.acos(u * 2 - 1) : Math.acos(Math.exp(Math.log(1 - u) * power) * 2 - 1)) / Math.PI;
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
    b2 = Math.tanh(Math.acosh(Math.cosh(a) / Math.cosh(b)) / 2),
    rot1 = 360 / p * DEGREE,
    rot2 = 360 / q * DEGREE,
    rot3 = 360 / r * DEGREE;
  let r01 = C(Math.cos(rot1), Math.sin(rot1)),
    r02 = C(Math.cos(rot2), Math.sin(rot2)),
    r03 = C(Math.cos(rot3), Math.sin(rot3));

  const c1 = C(h, 0),
    c2 = C(0, b1),
    c3 = C(0, b2);

  return z => {
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
  let result = 0;

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
    let nx1 = Math.cos(b * z.im);
    let nx2 = Math.sin(c * z.im);
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
    const xoff = z.re > 0 ? real : -real,
      yoff = z.im > 0 ? imaginary : -imaginary;
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

/* 3D */
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
    n1[0] + n2[0] + n3[0] + n4[0],
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
  let cosa = Math.cos(yaw);
  let sina = Math.sin(yaw);

  let cosb = Math.cos(pitch);
  let sinb = Math.sin(pitch);

  let cosc = Math.cos(roll);
  let sinc = Math.sin(roll);

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

function perspective3D(n, f) {
  let perspectiveMatrix = matrixMultiply(IDENTITY_MATRIX, matrixPerspectiveProjection(n, f));
  return z => {
    let result = applyMatrix([z.re, z.im, z.z], perspectiveMatrix);
    return {
      ...z,
      re: result[0],
      im: result[1],
      z: -Math.sqrt(Math.hypot(z.re, z.im, z.z))
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

function heightMap(zRange) {
  return z => {
    let height = z.z === 0 ? 0 : (z.z + zRange) / (2 * zRange);
    return {
      ...z,
      red: height,
      green: height,
      blue: height,
    }
  }
}

function basicLighting(theta, diffuse) {
  let sinTheta = Math.sin(theta * DEGREE);
  let cosTheta = Math.cos(theta * DEGREE);

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
    let brightness = z.z === 0 ? 0 : halfLength / (halfLength + Math.max(0, (z.z - startZ) / startZ));
    return {
      ...z,
      red: z.red * brightness + (1 - brightness) * mistColor.red,
      green: z.green * brightness + (1 - brightness) * mistColor.green,
      blue: z.blue * brightness + (1 - brightness) * mistColor.blue,
    }
  }
}

// TODO: maybe use textures instead
function enviormentLight(x, y, z, lights) {
  let left = 1;
  let sum = 0;
  for(let i = 0; i < lights.length; i++) {
    let d = Math.hypot(x - lights[i][0][0], y + lights[i][0][1], z - lights[i][0][2]);
    if(d < lights[i][1] + Math.SQRT2) {
      let t = Math.min(1, (lights[i][1] + Math.SQRT2 - d) / Math.SQRT2);
      sum += left * t * lights[i][2];;
      left -= left * t;
    }
  }
  return sum;
}

function lightRoomEnviorment(x, y, z) {
  let lights = [
    [normalize3([3, 3, 1]), 0.4, 1],
    [normalize3([-3, 3, 1]), 0.4, 1],
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
  let backdrop = normalize3([0, -1, -1]);
  if(y > 0.2 || Math.hypot(x - backdrop[0], y + backdrop[1], z - backdrop[2]) < 1.2) {
    let c = Math.max(0.2, y) * 0.5 + 0.3;
    return {
      red: c,
      green: c,
      blue: c
    };
  }
  if(y < -0.3) {
    return {
      red: y * 0.1 + 0.2,
      green: y * 0.1 + 0.2,
      blue: y * 0.1 + 0.2
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
  return enviormentLight(x, y, z, [
    [normalize3([3, 3, 1]), 0.4, 1],
    [normalize3([-3, 3, 1]), 0.4, 1],
  ]);
}

function dayEnviorment(x, y, z) {
  let lights = [
    [normalize3([-1, 2, 1]), 0.2, 1],
    [normalize3([3, 2, 2]), 0.2, 0.9],
    [normalize3([3.75, 2.5, 0.5]), 0.3, 0.9],
    [normalize3([3.5, 2, -0.5]), 0.23, 0.9],
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
  if(y < 0.2) {
    let c = -y * 0.2 + 0.8;
    return {
      red: c * 0.8,
      green: c * 0.9,
      blue: c * 1
    }
  }
  if(y < 0.3) {
    let c = -y * 0.2 + 0.7;
    return {
      red: c * 0.8,
      green: c * 0.8,
      blue: c * 1
    }
  }
  let c1 = y;
  let c2 = 1 - y;
  return {
    red: c1 * 0.3 + c2 * 0.5,
    green: c1 * 0.7 + c2 * 0.4,
    blue: c1 * 0.2 + c2 * 0.7
  };
}

function dayLights(x, y, z) {
  return enviormentLight(x, y, z, [
    [normalize3([-1, 2, 1]), 0.2, 1],
    [normalize3([3, 2, 2]), 0.2, 0.4],
    [normalize3([3.75, 2.5, 0.5]), 0.3, 0.4],
    [normalize3([3.5, 2, -0.5]), 0.23, 0.4],
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

function advancedLighting(theta1, theta2, ior, enviorment) {
  let skyBox = [lightRoomEnviorment, dayEnviorment][enviorment];
  let skyBoxLights = [lightRoomLights, dayLights][enviorment];
  let rotation1 = rotate3D(theta1, 0, 0);
  let rotation2 = rotate3D(0, theta2, 0);

  return z => {
    if(z.z === 0){
      return z;
    }
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
        ...lerp(
          {red: 0, green: 0, blue: 0},
          lerp(
            z,
            skyBox(n.re, n.im, n.z),
            schlick(ior, _normal)),
          brightness)
      };
    }
    let amp = Math.sqrt(normal[0] * normal[0] + normal[1] * normal[1]);
    normal = [normal[0] / amp, normal[1] / amp, (normal[2] / amp) / Math.min(z.width, z.height) * Math.abs(z.z)];
    let _at = [z.re, z.im, z.z];

    while(true) {
      _at = [_at[0] + normal[0], _at[1] + normal[1], _at[2] + normal[2]];
      let at = [Math.round(_at[0]), Math.round(_at[1]), _at[2]];
      if(at[0] < 0 || at[1] < 0 || at[0] + 1 >= z.width || at[1] + 1 >= z.height) {
        return {
          ...z,
          ...lerp(
            {red: 0, green: 0, blue: 0},
            lerp(
              z,
              skyBox(n.re, n.im, n.z),
              schlick(ior, _normal)),
            brightness)
        };
      }
      let sample = z.zBuffer[at[0] + at[1] * z.width];
      if(sample !== 0 && sample > at[2]) {
        let bounceNormal = getNormal({
          ...z,
          re: at[0],
          im: at[1],
          z: sample
        })
        newN = rotation1(rotation2({
          re: bounceNormal[0],
          im: bounceNormal[1],
          z: bounceNormal[2]
        }));
        let bounceBrightness = skyBoxLights(n.re, n.im, n.z);
        return {
          ...z,
          ...lerp(
            {red: 0, green: 0, blue: 0},
            lerp(
              z,
              lerp(
                {red: 0, green: 0, blue: 0},
                lerp(
                  {
                    red: z.mainBuffer[(at[0] + at[1] * z.width) * 3] / 255,
                    green: z.mainBuffer[(at[0] + at[1] * z.width) * 3 + 1] / 255,
                    blue: z.mainBuffer[(at[0] + at[1] * z.width) * 3 + 2] / 255
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

/* description s*/
const BUILT_IN_TRANSFORMS = {
  //shaders
  gamma: gamma,
  normalMap: normalMap,
  heightMap: heightMap,
  basicLighting: basicLighting,
  advancedLighting: advancedLighting,
  mist: mist,
  //3D transforms
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
  disc: disc,
  dragon: dragon,
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
  murl2: murl2,
  nSplit: nSplit,
  pdj: pdj,
  pointSymmetry: pointSymmetry,
  rotate: rotate,
  scale: scale,
  scale2: scale2,
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
  //shaders
  gamma: [{
    name: "gamma",
    type: "number",
    default: 2.2
  }],
  normalMap: [],
  heightMap: [{
    name: "z range",
    type: "number",
    default: 1
  }],
  basicLighting: [{
    name: "theta",
    type: "number",
    default: "90*DEGREE",
  }, {
    name: "diffuse",
    type: "number",
    default: 0.3,
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
    name: "enviorment",
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
  matrix3D: [{
    name: "matrix",
    type: "array",
    default: IDENTITY_MATRIX
  }],
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
  perspective3D: [{
    name: "n",
    type: "number",
    default: 1
  }, {
    name: "f",
    type: "number",
    default: 0
  }],
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
