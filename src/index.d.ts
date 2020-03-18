type Coords = { x: number; y: number; z?: number; } | { lat: number; lng: number; };
export type Beacon = { distance: number; } & Coords;

declare function locate(
  beacons: Beacon[],
  options?: {
    geometry?: '2d' | '3d' | 'earth',
    method?: 'lse' | 'lseInside' | 'circleSizing',
    maxIterations?: number
  }
): Coords;

export default locate;
