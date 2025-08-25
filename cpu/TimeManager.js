// --- Point class: Linked list node for trajectory points ---
class Point {
    static ID = 1;

    constructor(time, pt_data) {
        this.id = Point.ID++;
        this.time = time;
        this.data = pt_data;
        this.next = null;
        this.prev = null;
        this.observers = new Map(); // observerId -> this
    }

    CompareTime(otherTime) {
        return this.time - otherTime;
    }
}

// --- Trajectory: linked list of points with per-observer tracking ---
class Trajectory {
    constructor(id) {
        this.id = id;
        this.head = null;
        this.tail = null;
        this.length = 0;
        this.observers = new Map(); // observerId -> Point
    }

    Add(time, pt) {
        const node = new Point(time, pt);
        if (!this.head) {
            this.head = this.tail = node;
        } else {
            this.tail.next = node;
            node.prev = this.tail;
            this.tail = node;
        }
        this.length++;
    }

    Observed(observerId, node) {
        if (this.observers.has(observerId)) {
            const prev = this.observers.get(observerId);
            prev.observers.delete(observerId);
        }

        node.observers.set(observerId, node);
        this.observers.set(observerId, node);

        while (this.head && this.head.observers.size === 0) {
            const next = this.head.next;
            this.head.next = this.head.prev = null;
            this.head.observers.clear();
            this.head = next;
            if (next) next.prev = null;
            this.length--;
        }

        if (!this.head) this.tail = null;
    }

    StartNode(observerId) {
        return this.observers.get(observerId) || this.head;
    }

    IterateFrom(startNode, callback) {
        let node = startNode;
        while (node) {
            if (callback(node)) break;
            node = node.next;
        }
    }

    ForEach(callback) {
        this.IterateFrom(this.head, callback);
    }

    GetLatest() {
        return this.tail;
    }
}

// --- Trajectories: manages all particle trajectories ---
class Trajectories {
    constructor() {
        this.map = new Map(); // particleId -> Trajectory
    }

    Add(id, time, pt) {
        if (!this.map.has(id)) {
            const traj = new Trajectory(id);
            traj.Add(time, pt);
            this.map.set(id, traj);
        } else {
            this.map.get(id).Add(time, pt);
        }
    }

    Delete(id){
        this.map.delete(id);
    }

    Get(id) {
        return this.map.get(id);
    }

    Clear() {
        this.map.clear();
    }
}

// --- TimeManager: stores and queries observer-relative history ---
class TimeManager {
    constructor(max = 10000) {
        this.time_per_frame = Simulation.DistanceLightCoversInAFrame / PhysicsEngine.C;
        this.distance_per_frame = Simulation.DistanceLightCoversInAFrame;
        this.trajectories = new Trajectories();
    }

    Init(particles) {
        for (let observer of particles) {
            for (let observed of particles) {
                if (observer.id === observed.id) continue; // skip self
    
                // Ensure trajectory exists
                let traj = this.trajectories.Get(observed.id);
                if (!traj) {
                    traj = new Trajectory(observed.id);
                    this.trajectories.map.set(observed.id, traj);
                }
    
                // Initialize observer's start node to the head (null-safe)
                if (traj.head) {
                    traj.observers.set(observer.id, traj.head);
                    traj.head.observers.set(observer.id, traj.head);
                }
            }
        }
    }

    Store(particles, time) {
        for (const particle of particles) {
            this.trajectories.Add(particle.id, time, particle.frame_data);
        }
    }

    Delete(particles){
        for(const particle of particles){
            this.trajectories.Delete(particle.id);
        }
    }

    Get(observer, observed) {
        const observer_fd = observer.frame_data;
        const traj = this.trajectories.Get(observed.id);
        if (!traj || !traj.tail) return null;

        const observer_time = traj.tail.time;
        const start_node = traj.StartNode(observer.id) || traj.head;

        let result = null;

        traj.IterateFrom(start_node, point => {
            const observed_fd = point.data;

            const dx = observer_fd.x - observed_fd.x;
            const dy = observer_fd.y - observed_fd.y;
            const r_sqr = dx ** 2 + dy ** 2;
            const r = Math.sqrt(r_sqr);
            const delta = Math.abs((observer_time - point.time) - r / PhysicsEngine.C);

            if (delta < this.time_per_frame) {
                traj.Observed(observer.id, point);
                result = {
                    dx,
                    dy,
                    r,
                    r_sqr,
                    unit_x: dx / r,
                    unit_y: dy / r,
                    observed_data: observed_fd
                };
                return true; // stop iteration
            }
            return false;
        });

        return result;
    }

    Clear() {
        this.trajectories.Clear();
    }

    DrawPositionHistories(ctx, transformer) {
        ctx.save();
        ctx.lineWidth = 2;

        for (const [id, traj] of this.trajectories.map.entries()) {
            let first = true;
            ctx.beginPath();

            traj.ForEach(point => {
                const { x, y, color, type } = point.data;
                if (type !== "Photon") {
                    const [tx, ty] = transformer(x, y);

                    if (first) {
                        ctx.moveTo(tx, ty);
                        ctx.strokeStyle = color || "#ffffff";
                        first = false;
                    } else {
                        ctx.lineTo(tx, ty);
                    }
                }
            });

            if (!first) ctx.stroke();
        }

        ctx.restore();
    }

    DrawForcesHistory(ctx, transformer, maxFrames = 1000) {
        const ringSpacing = Simulation.NumPixelsLightCoversInAFrame * Simulation.ParticleForceRingSpacing;
        ctx.save();

        for (const [id, traj] of this.trajectories.map.entries()) {
            let node = traj.tail;
            let count = 0;

            while (node && count < maxFrames) {
                if (count % Simulation.ParticleForceRingSpacing === 0) {
                    const { x, y, q, type } = node.data;
                    if (type !== "Photon") {
                        const [px, py, pr] = transformer(x, y, count * Simulation.DistanceLightCoversInAFrame);

                        ctx.beginPath();
                        ctx.arc(px, py, pr, 0, Math.PI * 2, false);
                        ctx.lineWidth = ringSpacing+2;

                        const a = ((maxFrames - count) / maxFrames) ** 15;
                        ctx.strokeStyle = `rgba(${q >= 0 ? 255 : 0},${q <= 0 ? 255 : 0},255,${a})`;
                        ctx.stroke();
                    }
                }

                node = node.prev;
                count++;
            }
        }

        ctx.restore();
    }
}
