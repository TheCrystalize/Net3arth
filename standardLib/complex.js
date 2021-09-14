function C(real, imaginary){
  return {re:real, im: imaginary};
}

function compAdd(z, c) {
  return {
    re: z.re + c.re,
    im: z.im + c.im
  }
}

function compAddScalar(z, s) {
  return {
    re: z.re + s,
    im: z.im
  }
}

function compSub(z, c) {
  return {
    re: z.re - c.re,
    im: z.im - c.im
  }
}

function compMult(z, c) {
  return {
    re: z.re * c.re - z.im * c.im,
    im: z.re * c.im + z.im * c.re
  }
}

function compMultScalar(z, s) {
  return {
    re: z.re * s,
    im: z.im * s
  }
}

function compDiv(z, c) {
  const s = 1 / (c.re * c.re + c.im * c.im);
  return {
    re: (z.re * c.re + z.im * c.im) * s,
    im: (z.im * c.re - z.re * c.im) * s
  }
}

function compPow(z, p) {
  const n = p * Math.atan2(z.im, z.re)
  return {
    re: Math.cos(n),
    im: Math.sin(n) * Math.exp(p * Math.log(z.re * z.re + z.im * z.im))
  }
}

function compLog(z) {
  return {
    re: 0.5 * Math.log(z.re * z.re + z.im * z.im),
    im: Math.atan2(z.im, z.re)
  }
}

function compExp(z) {
  const e = exp(z.re)
  return {
    re: Math.cos(z.im) * e,
    im: Math.sin(z.im) * e
  }
}
