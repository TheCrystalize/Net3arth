/*helper functions*/
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
    re: Math.cos(n),
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

function sinh(z) {
  return add(divScalar(exp(z), 2), neg(divScalar(exp(neg(z)), 2)))
}

function cosh(z) {
  return add(divScalar(exp(neg(z)), 2), divScalar(exp(z), 2))
}

function tanh(z) {
  return div(sinh(z), cosh(z))
}
