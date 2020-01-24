import { minimize_Powell } from 'optimization-js';
import minimize_cobyla from './3rdparty/jscobyla/index';

import { Earth, Point, Circle } from './geometry';

const sum = vals => vals.reduce((agg, v) => agg + v, 0);

const dist = (a, b, mode='2d') => {
  if (mode === '2d') return ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2) ** .5;
  if (mode === '3d') return ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2) ** .5;
  if (mode === 'earth') return Earth.gcd(a[0], a[1], b[0], b[1]);
};

const sumErrors = (x, pnts, radii, mode) => {
  let e = 0;
  for (let i = 0, l = pnts.length; i < l; i++) {
    e += (dist(x, pnts[i].std(), mode) - radii[i]) ** 2;
  }
  return e;
};

const isDisjoint = (circles, isOnEarthSurface) => {
  // does not support sophisticated checking for disjoint area on earth surface mode
  for (let i = 0, l = circles.length; i < l; i++) {
    for (let j = i+1; j < l; j++) {
      if (!circles[j].touch(circles[i], isOnEarthSurface)) return true;
    }
  }
  return false;
};

function lseFind(circles, mode='2d', constrain = false) {
  const numPnts = circles.length;
  const radii = circles.map(p => p.r);
  const pnts = circles.map(p => p.c);
  const sumR = sum(radii);
  const weights = radii.map(r => (sumR - r) / ((numPnts - 1) * sumR));

  let p0 = new Point(0, 0, 0); // Starting point
  for (let i = 0; i < numPnts; i++) {
    p0 = p0.sum(weights[i]).sum(pnts[i]);
  }

  const x0 = p0.std();

  let answer;
  if (constrain) {
    // console.info('GC-LSE geolocating...');
    if (!isDisjoint(circles, mode === 'earth')) {
      const constrainFn = (x, beaconIdx) => radii[beaconIdx] - dist(x, pnts[beaconIdx].std(), mode);

      let x = x0.slice(); // mutated result
      minimize_cobyla(
        (n, m, x, con) => { // objective function
          for (let i = 0; i < m; i++) {
            con[i] = constrainFn(x, i); // one constrain per beacon
          }
          return sumErrors(x, pnts, radii, mode);
        },
        x.length, // # vars
        numPnts, // # constraints
        x0.slice(), // result
        1, // rhobeg
        1e-5, // rhoend
        0, // iprint
        1000 // maxfun
      );

      answer = x;
    } else {
      throw new Error('Disjointed beacons', circles);
    }
  } else {
    // console.info('LSE geolocating...');
    const res = minimize_Powell(x => sumErrors(x, pnts, radii, mode), x0);
    answer = res.argument;
  }

  return new Point(...answer);
}

function locate(beacons, { mode = '2d', constrain = false } = {}) {
  if (mode !== '2d' && mode !== '3d' && mode !== 'earth') {
    throw new Error(`Mode not supported: ${mode}`);
  }

  const circles = beacons.map(({ distance, ...coords }) => {
    const c = mode === '2d'
      ? new Point(coords.x, coords.y)
      : mode === '3d'
        ? new Point(coords.x, coords.y, coords.z)
        : new Point(coords.lng, coords.lat);

    return new Circle(c, distance); // in Earth mode, gcd distances are specified in meters
  });

  const { x, y, z } = lseFind(circles, mode, constrain);

  if (mode === 'earth') {
    let lng = x;
    let lat = y;

    // make lat, lngs uniform
    while (lng > 180) { lng -= 360 }
    while (lng < -180) { lng += 360 }
    while (lat > 90) { lat -= 180 }
    while (lat < -90) { lat += 180 }

    return { lat, lng };
  }

  return mode === '2d'
    ? { x, y }
    : { x, y, z };
}

export default locate;
