const std = @import("std");
const zigimg = @import("zigimg");

pub fn Vec(comptime size: usize, comptime T: type) type {
    return struct {
        const Self = @This();
        const Data = @Vector(size, T);

        data: Data,

        pub fn x(self: Self) T {
            return self.data[0];
        }
        pub fn y(self: Self) T {
            return self.data[1];
        }
        pub fn z(self: Self) T {
            return self.data[2];
        }

        pub fn maxDim(self: Self) u8 {
            var max_idx: u8 = 0;
            inline for (1..size) |i| {
                if (self.data[max_idx] < self.data[i]) {
                    max_idx = i;
                }
            }
            return max_idx;
        }

        pub fn zeroes() Self {
            return .{ .data = .{0, 0, 0}};
        }

        pub fn ones() Self {
            return .{ .data = .{1, 1, 1}};
        }

        pub fn init(_x: T, _y: T, _z: T) Self {
            return .{ .data = .{_x, _y, _z}};
        }

        pub fn fromArray(array: [size]T) Self {
            return .{ .data = array};
        }

        pub fn fromScalar(s: T) Self {
            return .{ .data = .{s, s, s}};
        }

        pub fn max(a: Self, b: Self) Self {
            return .{ .data = @max(a.data, b.data) };
        }

        pub fn min(a: Self, b: Self) Self {
            return .{ .data = @min(a.data, b.data) };
        }

        pub fn clamp(self: Self, min_value: T, max_value: T) Self {
            return min(self, max(fromScalar(min_value), fromScalar(max_value)));
        }

        pub fn scale(self: Self, s: T) Self {
            return .{ .data = self.data * @as(Data, @splat(s)) };
        }

        pub usingnamespace switch (@typeInfo(T)) {
            .Int => struct {
                pub fn toFloat(self: Self, comptime U: type) Vec(size, U) {
                    // TODO: https://github.com/ziglang/zig/issues/16267
                    return .{
                        .data = .{
                            @as(U, @floatFromInt(self.data[0])),
                            @as(U, @floatFromInt(self.data[1])),
                            @as(U, @floatFromInt(self.data[2])),
                        }
                    };
                }
            },
            .Float => struct {
                pub fn toInt(self: Self, comptime U: type) Vec(size, U) {
                    // TODO: https://github.com/ziglang/zig/issues/16267
                    return .{
                        .data = .{
                            @as(U, @intFromFloat(self.data[0])),
                            @as(U, @intFromFloat(self.data[1])),
                            @as(U, @intFromFloat(self.data[2])),
                        }
                    };
                }

                pub fn posInf() Self {
                    return .{
                        .data = .{
                            std.math.inf(T),
                            std.math.inf(T),
                            std.math.inf(T),
                        }
                    };
                }

                pub fn negInf() Self {
                    return .{
                        .data = .{
                            -std.math.inf(T),
                            -std.math.inf(T),
                            -std.math.inf(T),
                        }
                    };
                }

                pub fn length(self: Self) T {
                    return @sqrt(@reduce(.Add, self.data * self.data));
                }

                pub fn normalize(self: Self) Self {
                    return self.scale(1.0 / self.length());
                }

                pub fn abs(self: Self) Self {
                    return .{.data = @fabs(self.data)};
                }

                pub fn ceil(self: Self) Self {
                    return .{ .data = @ceil(self.data) };
                }

                threadlocal var prng = std.rand.DefaultPrng.init(0);

                pub fn randomUnitVector() Self {
                    // Using Gaussian distribution for all three coordinates of the vector
                    // will ensure an uniform distribution on the surface of the sphere.
                    const random = prng.random();
                    return init(
                        random.floatNorm(T),
                        random.floatNorm(T),
                        random.floatNorm(T),
                    ).normalize();
                }

                pub usingnamespace if (size == 3) struct {
                    pub fn toRGB(self: Self) zigimg.color.Rgb24 {
                        const rgb = self.clamp(0.0, 0.999999).scale(256);
                        return .{
                            .r = @intFromFloat(rgb.data[0]),
                            .g = @intFromFloat(rgb.data[1]),
                            .b = @intFromFloat(rgb.data[2]),
                        };
                    }
                } else struct {};
            },
            .Bool => struct {
                pub fn select(self: Self, a: anytype, b: anytype) @TypeOf(a, b) {
                    const Result = @TypeOf(a, b);
                    const Elem = std.meta.Child(std.meta.FieldType(Result, .data));
                    return .{.data = @select(Elem, self.data, a.data, b.data)};
                }
            },
            else => struct {}
        };

        pub usingnamespace if (size == 3) struct {
            pub fn cross(a: Self, b: Self) Self {
                const tmp0 = @shuffle(f32, a.data, a.data ,@Vector(3, i32){1,2,0});
                const tmp1 = @shuffle(f32, b.data, b.data ,@Vector(3, i32){2,0,1});
                const tmp2 = @shuffle(f32, a.data, a.data ,@Vector(3, i32){2,0,1});
                const tmp3 = @shuffle(f32, b.data, b.data ,@Vector(3, i32){1,2,0});
                return .{ .data = tmp0*tmp1-tmp2*tmp3 };
            }
        } else struct {};

        pub fn add(self: Self, b: Self) Self {
            return .{.data = self.data + b.data};
        }

        pub fn subtract(self: Self, b: Self) Self {
            return .{.data = self.data - b.data};
        }

        pub fn dot(a: Self, b: Self) T {
            return @reduce(.Add, a.data * b.data);
        }

        pub fn reduceMul(self: Self) T {
            return @reduce(.Mul, self.data);
        }

        pub fn mul(a: Self, b: Self) Self {
            return .{ .data = a.data * b.data };
        }

        pub fn div(a: Self, b: Self) Self {
            return .{ .data = a.data / b.data };
        }

        pub fn inc(self: Self) Self {
            return self.add(Self.fromScalar(1));
        }

        pub fn dec(self: Self) Self {
            return self.subtract(Self.fromScalar(1));
        }

        pub fn lessThan(a: Self, b: Self) Vec(size, bool) {
            return .{.data = a.data < b.data};
        }

        pub fn jsonParse(allocator: std.mem.Allocator, source: anytype, options: std.json.ParseOptions) !Self {
            return .{.data = try std.json.innerParse(Data, allocator, source, options)};
        }
    };
}


pub const Vec3 = Vec(3, f32);
pub const Vec3u = Vec(3, u32);

const vec3 = Vec3.init;
const vec3u = Vec3u.init;

test "cross product" {
    const a = vec3(1,-8,12);
    const b = vec3(4,6,3);
    const result = vec3(-96,45,38);
    try std.testing.expectEqual(Vec3.cross(a,b), result);
}

test "vector length" {
    const v = vec3(1.5, 100.0, -21.1);
    try std.testing.expectApproxEqAbs(v.length(), 102.21281720019266, 0.0001);
}

pub const Mat4 = struct {
    data: [4][4]f32,

    fn get(self: Mat4, row: usize, column: usize) f32 {
        return self.data[column][row];
    }

    fn set(self: *Mat4, row: usize, column: usize, val: f32) void {
        self.data[column][row] = val;
    }

    pub fn col3(self: Mat4, column: usize) Vec3 {
        return vec3(
            self.get(0, column),
            self.get(1, column),
            self.get(2, column),
        );
    }

    pub fn transformPosition(self: Mat4, v: Vec3) Vec3 {
        return .{ .data =
            self.col3(0).scale(v.x()).data +
            self.col3(1).scale(v.y()).data +
            self.col3(2).scale(v.z()).data +
            self.col3(3).data
        };
    }

    pub fn transformDirection(self: Mat4, v: Vec3) Vec3 {
        return .{ .data =
            self.col3(0).scale(v.x()).data +
            self.col3(1).scale(v.y()).data +
            self.col3(2).scale(v.z()).data
        };
    }
};

pub const Ray = struct {
    orig: Vec3,
    dir: Vec3,

    pub fn at(self: Ray, t: f32) Vec3 {
        return self.orig.add(self.dir.scale(t));
    }
};

const dot = Vec3.dot;
const cross = Vec3.cross;
const subtract = Vec3.subtract;
const add = Vec3.add;

pub const Bbox = struct {
    min: Vec3 = Vec3.posInf(),
    max: Vec3 = Vec3.negInf(),

    pub fn extendBy(self: Bbox, pos: Vec3) Bbox {
        return .{
            .min = Vec3.min(self.min, pos),
            .max = Vec3.max(self.max, pos),
        };
    }

    pub fn unionWith(a: Bbox, b: Bbox) Bbox {
        return .{
            .min = Vec3.min(a.min, b.min),
            .max = Vec3.max(a.max, b.max),
        };
    }

    pub fn size(self: Bbox) Vec3 {
        return subtract(self.max, self.min);
    }

    pub fn center(self: Bbox) Vec3 {
        return add(self.max, self.min).scale(0.5);
    }

    pub fn rayIntersection(self: Bbox, ray: Ray, t: *f32) bool {
        const sign = ray.dir.lessThan(Vec3.zeroes());

        const min = sign.select(self.max, self.min).subtract(ray.orig).div(ray.dir);
        const max = sign.select(self.min, self.max).subtract(ray.orig).div(ray.dir);

        var tmin = min.x();
        var tmax = max.x();

        if ((tmin > max.y()) or (tmax < min.y()))
            return false;

        tmin = @max(tmin, min.y());
        tmax = @min(tmax, max.y());

        if ((tmin > max.z()) or (tmax < min.z()))
            return false;

        tmin = @max(tmin, min.z());
        tmax = @min(tmax, max.z());

        t.* = tmin;

        return true;
    }
};

test "bbox rayIntersection 1" {
    const bbox = Bbox{
        .min = Vec3.fromScalar(-1),
        .max = Vec3.fromScalar(1),
    };
    const ray = Ray{
        .orig = vec3(0, 0, 5),
        .dir = vec3(0, 0, -1),
    };
    var t: f32 = undefined;
    try std.testing.expect(bbox.rayIntersection(ray, &t));
    try std.testing.expectApproxEqAbs(@as(f32, 4), t, 0.0001);
}

test "bbox rayIntersection 2" {
    const bbox = Bbox{
        .min = Vec3.fromScalar(1),
        .max = Vec3.fromScalar(2),
    };
    const ray = Ray{
        .orig = vec3(0, 0, 0),
        .dir = vec3(1, 1, 1).normalize(),
    };
    var t: f32 = undefined;
    try std.testing.expect(bbox.rayIntersection(ray, &t));
    try std.testing.expectApproxEqAbs(@sqrt(@as(f32, 3)), t, 0.0001);
}

test "bbox rayIntersection 3 (if orig is inside bbox => t < 0)" {
    const bbox = Bbox{
        .min = Vec3.fromScalar(-1),
        .max = Vec3.fromScalar(3),
    };
    const ray = Ray{
        .orig = vec3(0, 0, 0),
        .dir = vec3(1, 1, 0).normalize(),
    };
    var t: f32 = undefined;
    try std.testing.expect(bbox.rayIntersection(ray, &t));
    try std.testing.expectApproxEqAbs(-@sqrt(@as(f32, 2)), t, 0.0001);
}

test "bbox rayIntersection 4 (miss)" {
    const bbox = Bbox{
        .min = Vec3.fromScalar(-1),
        .max = Vec3.fromScalar(3),
    };
    const ray = Ray{
        .orig = vec3(5, 5, 5),
        .dir = vec3(1, 1, 0).normalize(),
    };
    var t: f32 = undefined;
    try std.testing.expect(bbox.rayIntersection(ray, &t) == false);
}

pub const Grid = struct {
    bbox: Bbox,
    resolution: Vec3u,
    cell_size: Vec3,

    pub fn init(bbox: Bbox, resolution: Vec3u) Grid {
        return .{
            .bbox = bbox,
            .resolution = resolution,
            .cell_size = Vec3.div(bbox.size(), resolution.toFloat(f32)),
        };
    }

    pub fn numCells(grid: Grid) u32 {
        return grid.resolution.reduceMul();
    }

    pub fn getCellPos(grid: Grid, point: Vec3) Vec3u {
        const pos = point.subtract(grid.bbox.min).div(grid.cell_size).toInt(u32);
        return Vec3u.min(pos, grid.resolution.dec());
    }

    pub fn getCellIdx(grid: Grid, x: usize, y: usize, z: usize) usize {
        return z * grid.resolution.x() * grid.resolution.y() + y * grid.resolution.x() + x;
    }

    pub fn traceRay(grid: Grid, ray: Ray) ?Iterator {
        var t_hit: f32 = undefined;
        if (grid.bbox.rayIntersection(ray, &t_hit) == false) {
            return null;
        }
        t_hit = @max(0, t_hit);

        const sign = ray.dir.lessThan(Vec3.zeroes());
        const step = sign.select(Vec3u.fromScalar(std.math.maxInt(u32)), Vec3u.fromScalar(1));
        const exit = sign.select(Vec3u.fromScalar(0), grid.resolution.dec());

        const t_delta = grid.cell_size.div(ray.dir).abs();

        const hit_local_pos = ray.at(t_hit).subtract(grid.bbox.min);
        const cell = Vec3u.min(hit_local_pos.div(grid.cell_size).toInt(u32), grid.resolution.dec());
        const next_cell = cell.add(sign.select(Vec3u.fromScalar(0), Vec3u.fromScalar(1))).toFloat(f32);
        const t_next_crossing = Vec3.fromScalar(t_hit).add(next_cell.mul(grid.cell_size).subtract(hit_local_pos).div(ray.dir));

        var it: Iterator = .{
            .cell = cell.data,
            .exit = exit.data,
            .step = step.data,
            .t_delta = t_delta.data,
            .t_next_crossing = t_next_crossing.data,
        };
        return it;
    }

    pub const Iterator = struct {
        cell: [3]u32,
        exit: [3]u32,
        step: [3]u32,
        t_delta: [3]f32,
        t_next_crossing: [3]f32,

        pub fn next(self: *Iterator) f32 {
            const k: u3 =
                (@as(u3, @intCast(@intFromBool(self.t_next_crossing[0] < self.t_next_crossing[1]))) << 2) +
                (@as(u3, @intCast(@intFromBool(self.t_next_crossing[0] < self.t_next_crossing[2]))) << 1) +
                (@as(u3, @intCast(@intFromBool(self.t_next_crossing[1] < self.t_next_crossing[2]))));
            const map = [8]u2{ 2, 1, 2, 1, 2, 2, 0, 0 };
            const axis = map[k];

            if (self.cell[axis] == self.exit[axis]) {
                return std.math.inf(f32);
            }

            const t_next_crossing = self.t_next_crossing[axis];

            self.cell[axis] = @addWithOverflow(self.cell[axis], self.step[axis])[0];
            self.t_next_crossing[axis] += self.t_delta[axis];

            return t_next_crossing;
        }
    };
};

test "decrement via add" {
    var x: u32 = 5;
    x = @addWithOverflow(x, std.math.maxInt(u32))[0];
    try std.testing.expectEqual(x, 4);
}

test "grid traceRay 1" {
    const grid = Grid.init(.{
        .min = vec3(0,0,0),
        .max = vec3(5,5,5),
    }, vec3u(5,5,5));
    const ray = .{
        .orig = vec3(0.5, 0.5, 0.5),
        .dir = vec3(2, 1, 0).normalize(),
    };
    var it = grid.traceRay(ray).?;
    try std.testing.expectEqual(it.cell, .{0, 0, 0});
    try std.testing.expectApproxEqAbs(it.next(), 0.559017002, 0.0001);
    try std.testing.expectEqual(it.cell, .{1, 0, 0});
    try std.testing.expectApproxEqAbs(it.next(), 1.11803400, 0.0001);
    try std.testing.expectEqual(it.cell, .{1, 1, 0});
    try std.testing.expectApproxEqAbs(it.next(), 1.67705106, 0.0001);
    try std.testing.expectEqual(it.cell, .{2, 1, 0});
    try std.testing.expectApproxEqAbs(it.next(), 2.79508495, 0.0001);
    try std.testing.expectEqual(it.cell, .{3, 1, 0});
    try std.testing.expectApproxEqAbs(it.next(), 3.35410213, 0.0001);
    try std.testing.expectEqual(it.cell, .{3, 2, 0});
    try std.testing.expectApproxEqAbs(it.next(), 3.91311883, 0.0001);
    try std.testing.expectEqual(it.cell, .{4, 2, 0});
    try std.testing.expectEqual(it.next(), std.math.inf(f32));
}

test "grid traceRay 2" {
    const grid = Grid.init(.{
        .min = vec3(0,0,0),
        .max = vec3(5,5,5),
    }, vec3u(5,5,5));
    const ray = .{
        .orig = vec3(0.5, 10.0, 0.5),
        .dir = vec3(0,-1,0),
    };
    var it = grid.traceRay(ray).?;
    try std.testing.expectEqual(it.cell, .{0, 4, 0});
    try std.testing.expectApproxEqAbs(it.next(), 6, 0.0001);
    try std.testing.expectEqual(it.cell, .{0, 3, 0});
    try std.testing.expectApproxEqAbs(it.next(), 7, 0.0001);
    try std.testing.expectEqual(it.cell, .{0, 2, 0});
    try std.testing.expectApproxEqAbs(it.next(), 8, 0.0001);
    try std.testing.expectEqual(it.cell, .{0, 1, 0});
    try std.testing.expectApproxEqAbs(it.next(), 9, 0.0001);
    try std.testing.expectEqual(it.cell, .{0, 0, 0});
    try std.testing.expectEqual(it.next(), std.math.inf(f32));
}

test "grid traceRay 3" {
    const grid = Grid.init(.{
        .min = vec3(0,0,0),
        .max = vec3(5,5,5),
    }, vec3u(5,5,5));
    const ray = .{
        .orig = vec3(0.5, -5.0, 0.5),
        .dir = vec3(0,1,0),
    };
    var it = grid.traceRay(ray).?;
    try std.testing.expectEqual(it.cell, .{0, 0, 0});
    try std.testing.expectApproxEqAbs(it.next(), 6, 0.0001);
    try std.testing.expectEqual(it.cell, .{0, 1, 0});
    try std.testing.expectApproxEqAbs(it.next(), 7, 0.0001);
    try std.testing.expectEqual(it.cell, .{0, 2, 0});
    try std.testing.expectApproxEqAbs(it.next(), 8, 0.0001);
    try std.testing.expectEqual(it.cell, .{0, 3, 0});
    try std.testing.expectApproxEqAbs(it.next(), 9, 0.0001);
    try std.testing.expectEqual(it.cell, .{0, 4, 0});
    try std.testing.expectEqual(it.next(), std.math.inf(f32));
}

test "grid traceRay 4" {
    const grid = Grid.init(.{
        .min = vec3(0,0,0),
        .max = vec3(5,5,5),
    }, vec3u(5,5,5));
    const ray = .{
        .orig = vec3(0.5, 0.5, 0.5),
        .dir = vec3(1, 1, 0).normalize(),
    };
    var it = grid.traceRay(ray).?;
    try std.testing.expectEqual(it.cell, .{0, 0, 0});
    try std.testing.expectApproxEqAbs(it.next(), 0.707106769, 0.0001);
    try std.testing.expectEqual(it.cell, .{0, 1, 0});
    try std.testing.expectApproxEqAbs(it.next(), 0.707106769, 0.0001);
    try std.testing.expectEqual(it.cell, .{1, 1, 0});
    try std.testing.expectApproxEqAbs(it.next(), 2.12132024, 0.0001);
    try std.testing.expectEqual(it.cell, .{1, 2, 0});
    try std.testing.expectApproxEqAbs(it.next(), 2.12132024, 0.0001);
    try std.testing.expectEqual(it.cell, .{2, 2, 0});
    try std.testing.expectApproxEqAbs(it.next(), 3.53553390, 0.0001);
    try std.testing.expectEqual(it.cell, .{2, 3, 0});
    try std.testing.expectApproxEqAbs(it.next(), 3.53553390, 0.0001);
    try std.testing.expectEqual(it.cell, .{3, 3, 0});
    try std.testing.expectApproxEqAbs(it.next(), 4.94974756, 0.0001);
    try std.testing.expectEqual(it.cell, .{3, 4, 0});
    try std.testing.expectApproxEqAbs(it.next(), 4.94974756, 0.0001);
    try std.testing.expectEqual(it.cell, .{4, 4, 0});
    try std.testing.expectEqual(it.next(), std.math.inf(f32));
}

pub const Triangle = struct {
    v0: Vec3,
    e1: Vec3,
    e2: Vec3,

    pub fn init(v0: Vec3, v1: Vec3, v2: Vec3) Triangle {
        return .{
            .v0 = v0,
            .e1 = subtract(v1, v0),
            .e2 = subtract(v2, v0),
        };
    }

    pub fn rayIntersection(self: Triangle, ray: Ray, t: *f32, _u: *f32, _v: *f32) bool
    {
        const pvec = cross(ray.dir, self.e2);
        const det = dot(self.e1, pvec);

        const epsilon = 0.00000001; // TODO

        // // if the determinant is negative, the triangle is 'back facing'
        // // if the determinant is close to 0, the ray misses the triangle
        if (det < epsilon) return false;

        const inv_det = 1.0 / det;

        const tvec = subtract(ray.orig, self.v0);
        const u = dot(tvec, pvec) * inv_det;
        if (u < 0 or u > 1) return false;

        const qvec = cross(tvec, self.e1);
        const v = dot(ray.dir, qvec) * inv_det;
        if (v < 0 or u+v > 1) return false;

        t.* = dot(self.e2, qvec) * inv_det;
        _u.* = u;
        _v.* = v;

        return true;
    }
};
