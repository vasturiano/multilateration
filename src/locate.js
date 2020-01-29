import { nelderMead as fmin } from 'fmin/index.js';
// import { intersectionArea } from 'venn.js';
import { intersectionArea, containedInCircles } from 'venn.js/src/circleintersection';

import { Earth, Point, Circle } from './geometry';

const METERS_PER_DEGREE = Earth.gcd(0, 0, 0, 1);

const sum = vals => vals.reduce((agg, v) => agg + v, 0);
const avg = vals => sum(vals) / vals.length;

const dist = (a, b, geometry = '2d') => {
  if (geometry === '2d')
    return ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2) ** 0.5;
  if (geometry === '3d')
    return (
      ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2) ** 0.5
    );
  if (geometry === 'earth') return Earth.gcd(a[0], a[1], b[0], b[1]);
};

const calcSumErrors = (coords, radiusProportion, circles, geometry) =>
  sum(
    circles.map(
      ({ c, r }) =>
        (dist(coords, c.std(), geometry) - r * radiusProportion) ** 2
    )
  );

const calcIntersectionArea = (radiusProportion, circles, geometry, stats) => {
  if (!['2d', 'earth'].includes(geometry)) {
    throw new Error(`Unsupported geometry: ${geometry}`);
  }

  const rFactor =
    radiusProportion / (geometry === 'earth' ? METERS_PER_DEGREE : 1);
  const circleData = circles.map(({ c, r }) => ({
    x: c.x,
    y: c.y,
    radius: r * rFactor
  }));

  return (
    intersectionArea(circleData, stats) *
    (geometry === 'earth' ? METERS_PER_DEGREE ** 2 : 1)
  );
};

const calcIntersectionCentroid = (radiusProportion, circles, geometry) => {
  const stats = {};
  calcIntersectionArea(radiusProportion, circles, geometry, stats);

  return ['x', 'y'].map(c => avg(stats.innerPoints.map(p => p[c])));
};

const calcStartingCoords = (circles, geometry) => {
  // center coordinates of smallest circle
  const smallestCircle = circles.slice().sort((a, b) => a.r - b.r)[0];
  const smallestCircleR =
    smallestCircle.r / (geometry === 'earth' ? METERS_PER_DEGREE : 1);
  const smallestCircleC = smallestCircle.c;

  // weighted centroid of all pnts
  const sumR = sum(circles.map(p => p.r));
  let wCentroid = new Point(0, 0, 0);
  circles.forEach(({ c, r }) => {
    const weight = (sumR - r) / ((circles.length - 1) * sumR);
    wCentroid = wCentroid.sum(c.mult(weight));
  });

  const p0 =
    smallestCircleC.dist(wCentroid) <= smallestCircleR
      ? wCentroid // pick weighted centroid if it's included within the smallest circle radius
      : smallestCircleC.sum(
          new Point(
            // jiggle it within the smallest circle radius
            ...[...new Array(3)].map(
              () => (Math.random() * 2 - 1) * smallestCircleR * 0.7
            )
          )
        );

  return p0.std().slice(0, geometry === '3d' ? 3 : 2);
};

const isDisjoint = (circles, isOnEarthSurface) => {
  // does not support sophisticated checking for disjoint area on earth surface mode
  for (let i = 0, l = circles.length; i < l; i++) {
    for (let j = i + 1; j < l; j++) {
      if (!circles[j].touch(circles[i], isOnEarthSurface)) return true;
    }
  }
  return false;
};

const isPntInOverlapArea = (pnt, circles, geometry) => {
  if (!['2d', 'earth'].includes(geometry)) {
    throw new Error(`Unsupported geometry: ${geometry}`);
  }

  const rFactor = 1 / (geometry === 'earth' ? METERS_PER_DEGREE : 1);

  const pntData = { x: pnt[0], y: pnt[1] };

  const circleData = circles.map(({ c, r }) => ({
    x: c.x,
    y: c.y,
    radius: r * rFactor
  }));

  return containedInCircles(pntData, circleData);
};

function findCoords(circles, { geometry = '2d', method = 'lse' }) {
  if (method === 'lse') {
    return fmin(
      coords => calcSumErrors(coords, 1, circles, geometry),
      calcStartingCoords(circles, geometry),
      { maxIterations: 100 }
    ).x;
  } else if (method === 'lseInside') {
    if (isDisjoint(circles, geometry === 'earth')) {
      throw new Error('Disjointed beacons', circles);
    }

    const startingCoords = calcStartingCoords(circles, geometry);

    const [, ...coordsLse] = fmin(
      ([radiusProportion, ...coords]) => {
        // don't enlarge circles to make it fit
        if (radiusProportion > 1) return Infinity;

        if (geometry === 'earth') {
          const [lng, lat] = coords;
          if (Math.abs(lng) > 180 || Math.abs(lat) > 90) return Infinity; // lat, lng out of bounds
        }

        if (!isPntInOverlapArea(coords, circles, geometry)) return Infinity;

        return calcSumErrors(coords, radiusProportion, circles, geometry);
      },
      [1, ...startingCoords],
      { maxIterations: 100 }
    ).x;

    return coordsLse;
  } else if (method === 'circleSizing') {
    const [minRadiusProportion] = fmin(
      ([radiusProportion]) => {
        const area = calcIntersectionArea(radiusProportion, circles, geometry);
        return area > 0 ? area : Infinity;
      },
      [1],
      { maxIterations: 100 }
    ).x;

    return calcIntersectionCentroid(minRadiusProportion, circles, geometry);
  }
}

function locate(beacons, { geometry = '2d', method = 'lse' } = {}) {
  if (!['2d', '3d', 'earth'].includes(geometry)) {
    throw new Error(`Unsupported geometry: ${geometry}`);
  }

  if (!['lse', 'lseInside', 'circleSizing'].includes(method)) {
    throw new Error(`Unsupported method: ${method}`);
  }

  const circles = beacons.map(({ distance, ...coords }) => {
    const c =
      geometry === '2d'
        ? new Point(coords.x, coords.y)
        : geometry === '3d'
        ? new Point(coords.x, coords.y, coords.z)
        : new Point(coords.lng, coords.lat);

    return new Circle(c, distance); // in Earth mode, gcd distances are specified in meters
  });

  const [x, y, z] = findCoords(circles, { geometry, method });

  if (geometry === 'earth') {
    let lng = x;
    let lat = y;

    // make lat, lngs uniform
    while (lng > 180) {
      lng -= 360;
    }
    while (lng < -180) {
      lng += 360;
    }
    while (lat > 90) {
      lat -= 180;
    }
    while (lat < -90) {
      lat += 180;
    }

    return { lat, lng };
  }

  return geometry === '2d' ? { x, y } : { x, y, z };
}

export default locate;
