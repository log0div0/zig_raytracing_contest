
const std = @import("std");
const Gltf = @import("zgltf");
const zigimg = @import("zigimg");

const main = @import("main.zig");
const linalg = @import("linalg.zig");
const stage2 = @import("stage2.zig");
const stage3 = @import("stage3.zig");

const Vec3 = linalg.Vec3;
const Mat4 = linalg.Mat4;
const Grid = linalg.Grid;
const Bbox = linalg.Bbox;

const vec3 = Vec3.init;
const cross = Vec3.cross;






// =====================================================================================
// =====================================================================================
// =====================================================================================
// =====================================================================================

fn loadGltfImagesWorker(tmp_allocator: std.mem.Allocator, arena_allocator: std.mem.Allocator,
    gltf_dir: []const u8, gltf: *Gltf,
    thread_idx: usize, thread_num: usize) !void
{
    var image_idx = thread_idx;
    while (image_idx < gltf.data.images.items.len) : (image_idx += thread_num) {
        const image = &gltf.data.images.items[image_idx];
        const im = blk: {
            if (image.buffer_view != null) {
                const buffer_view = gltf.data.buffer_views.items[image.buffer_view.?];
                const buffer = gltf.data.buffers.items[buffer_view.buffer];
                const begin = buffer_view.byte_offset;
                const end = begin + buffer_view.byte_length;
                break :blk try zigimg.Image.fromMemory(arena_allocator, buffer.data.?[begin..end]);
            } else {
                var tmp = [_]u8{undefined} ** 256;
                const img_path = try std.fmt.bufPrintZ(&tmp, "{s}/{s}", .{gltf_dir, image.uri.?});
                std.log.debug("Loading image {s}", .{img_path});
                const data = try main.loadFile(img_path, tmp_allocator);
                defer tmp_allocator.free(data);
                break :blk try zigimg.Image.fromMemory(arena_allocator, data);
            }
        };
        const ptr = try arena_allocator.create(zigimg.Image);
        ptr.* = im;
        image.data = @ptrCast(ptr);
    }
}

pub fn loadGltfFile(_allocator: std.mem.Allocator, path: []const u8, threads: []std.Thread) !Gltf {
    const gltf_dir = std.fs.path.dirname(path).?;

    var gltf = Gltf.init(_allocator);
    errdefer gltf.deinit();

    std.log.debug("Loading scene {s}", .{path});
    const arena_allocator = gltf.arena.allocator();
    const buf = try main.loadFile(path, arena_allocator);
    try gltf.parse(buf);

    for (gltf.data.buffers.items, 0..) |*buffer, i| {
        if (i == 0 and buffer.uri == null) {
            buffer.data = gltf.glb_binary;
            continue;
        }
        var tmp = [_]u8{undefined} ** 256;
        const buf_path = try std.fmt.bufPrint(&tmp, "{s}/{s}", .{gltf_dir, buffer.uri.?});
        std.log.debug("Loading buffer {s}", .{buf_path});
        buffer.data = try main.loadFile(buf_path, arena_allocator);
    }

    var safe_tmp_allocator = std.heap.ThreadSafeAllocator{.child_allocator = _allocator};
    var safe_arena_allocator = std.heap.ThreadSafeAllocator{.child_allocator = arena_allocator};
    for (threads, 0..) |*thread, i| {
        thread.* = try std.Thread.spawn(.{}, loadGltfImagesWorker, .{
            safe_tmp_allocator.allocator(), safe_arena_allocator.allocator(),
            gltf_dir, &gltf,
            i, threads.len
        });
    }
    for (threads) |*thread| {
        thread.join();
    }

    return gltf;
}






// =====================================================================================
// =====================================================================================
// =====================================================================================
// =====================================================================================

fn findPrimitiveAttribute(primitive: Gltf.Primitive, comptime tag: std.meta.Tag(Gltf.Attribute)) !?Gltf.Index {
    for (primitive.attributes.items) |attribute| {
        if (attribute == tag) {
            return @field(attribute, @tagName(tag));
        }
    }
    return error.AttributeNotFound;
}

fn Accessor(comptime T: type) type {
    return struct {
        const Self = @This();

        base_address: [*]const u8,
        num: usize,
        stride: usize,

        fn init(gltf: Gltf, accessor_idx: ?Gltf.Index) Self
        {
            const accessor = gltf.data.accessors.items[accessor_idx.?];
            switch (T) {
                Vec3 => {
                    std.debug.assert(accessor.type == .vec3);
                    std.debug.assert(accessor.component_type == .float);
                },
                [2]f32 => {
                    std.debug.assert(accessor.type == .vec2);
                    std.debug.assert(accessor.component_type == .float);
                },
                u16 => {
                    std.debug.assert(accessor.type == .scalar);
                    std.debug.assert(accessor.component_type == .unsigned_short);
                },
                else => {
                    @compileError("Implement me");
                }
            }

            const buffer_view = gltf.data.buffer_views.items[accessor.buffer_view.?];
            const buffer = gltf.data.buffers.items[buffer_view.buffer];

            return .{
                .base_address =
                    buffer.data.?.ptr +
                    buffer_view.byte_offset +
                    accessor.byte_offset,
                .num = accessor.count,
                .stride = accessor.stride,
            };
        }

        fn at(self: Self, idx: usize) T {
            switch (T) {
                Vec3 => {
                    const ptr: [*]const f32 = @ptrCast(@alignCast(self.base_address + idx*self.stride));
                    return Vec3.init(ptr[0],ptr[1],ptr[2]);
                },
                [2]f32 => {
                    const ptr: [*]const f32 = @ptrCast(@alignCast(self.base_address + idx*self.stride));
                    return ptr[0..2].*;
                },
                u16 => {
                    const ptr: [*]const u16 = @ptrCast(@alignCast(self.base_address + idx*self.stride));
                    return ptr[0];
                },
                else => {
                    @compileError("Implement me");
                }
            }
        }
    };
}

fn calcTriangles(gltf: Gltf) usize {
    var num_triangles: usize = 0;
    for (gltf.data.nodes.items) |node| {
        if (node.mesh != null) {
            const mesh = gltf.data.meshes.items[node.mesh.?];
            for (mesh.primitives.items) |primitive| {
                const accessor = gltf.data.accessors.items[primitive.indices.?];
                num_triangles += accessor.count / 3;
            }
        }
    }
    return num_triangles;
}

fn loadTriangles(gltf: Gltf, triangles: *std.MultiArrayList(stage2.Triangle)) !void {
    var counter: usize = 0;
    for (gltf.data.nodes.items) |node| {
        if (node.mesh != null) {
            const mesh = gltf.data.meshes.items[node.mesh.?];
            for (mesh.primitives.items) |primitive|
            {
                std.debug.assert(primitive.mode == .triangles);

                const positions = Accessor(Vec3).init(gltf, try findPrimitiveAttribute(primitive, .position));
                const normals = Accessor(Vec3).init(gltf, try findPrimitiveAttribute(primitive, .normal));
                const texcoords = Accessor([2]f32).init(gltf, try findPrimitiveAttribute(primitive, .texcoord));
                const indices = Accessor(u16).init(gltf, primitive.indices);

                const triangles_count = indices.num / 3;

                const matrix = Mat4{.data = gltf.getGlobalTransform(node)};

                for (0..triangles_count) |triangle_idx| {
                    var bbox: Bbox = .{};
                    var pos: [3]Vec3 = undefined;
                    var data = stage3.Triangle.Data {
                        .v = undefined,
                        .material_idx = primitive.material.?,
                    };
                    for (0..3) |i| {
                        const index_idx = triangle_idx*3+i;
                        const vertex_idx = indices.at(index_idx);
                        pos[i] = matrix.transformPosition(positions.at(vertex_idx));
                        bbox = bbox.extendBy(pos[i]);
                        data.v[i] = .{
                            .normal = matrix.transformDirection(normals.at(vertex_idx)).normalize(), // TODO: use adjusent matrix
                            .texcoord = texcoords.at(vertex_idx),
                        };
                    }
                    triangles.set(counter, .{
                        .centroid = bbox.center().data,
                        .bbox = bbox,
                        .pos = pos,
                        .data = data,
                    });
                    counter += 1;
                }
            }
        }
    }
}

pub fn loadGeometry(gltf: Gltf, geometry: *stage2.Geometry) !void {
    const allocator = geometry.arena.allocator();

    const num_triangles = calcTriangles(gltf);

    var triangles = std.MultiArrayList(stage2.Triangle){};
    try triangles.resize(allocator, num_triangles);

    try loadTriangles(gltf, &triangles);

    geometry.triangles = triangles.slice();
}




// =====================================================================================
// =====================================================================================
// =====================================================================================
// =====================================================================================

fn findCameraNode(gltf: Gltf) !Gltf.Node {
    for (gltf.data.nodes.items) |node| {
        if (node.camera != null) {
            return node;
        }
    }
    return error.CameraNodeNotFound; // TODO: implement recursive search / deal with multiple instances of the same camera
}

pub fn loadCamera(gltf: Gltf, width: ?u16, height: ?u16) !stage3.Camera {

    const camera_node = try findCameraNode(gltf);

    const camera = gltf.data.cameras.items[camera_node.camera.?];
    if (camera.type != .perspective) {
        return error.OnlyPerspectiveCamerasSupported;
    }

    var w: u32 = 0;
    var h: u32 = 0;

    if (width == null and height == null)
    {
        return error.OutputImgSizeIsNotSpecified;
    }
    else if (width != null and height != null)
    {
        if (camera.type.perspective.aspect_ratio != null) {
            return error.CameraHasAspectRatio;
        }
        w = width.?;
        h = height.?;
    }
    else
    {
        if (camera.type.perspective.aspect_ratio == null) {
            return error.CameraHasntAspectRatio;
        }
        const aspect_ratio = camera.type.perspective.aspect_ratio.?;
        w = width orelse @intFromFloat(@as(f32, @floatFromInt(height.?)) * aspect_ratio);
        h = height orelse @intFromFloat(@as(f32, @floatFromInt(width.?)) / aspect_ratio);
    }

    std.log.info("Pixels count: {}", .{w*h});

    const f_w: f32 = @floatFromInt(w);
    const f_h: f32 = @floatFromInt(h);

    const matrix = Mat4{.data = gltf.getGlobalTransform(camera_node)};
    const origin = matrix.col3(3);
    const fwd = matrix.col3(2).scale(-1).normalize();
    const world_up = vec3(0,1,0);

    const right = cross(fwd, world_up).normalize();
    const up = cross(fwd, right);

    const focal_length = (f_h / 2) / @tan(camera.type.perspective.yfov / 2);

    const lower_left_corner = fwd.scale(focal_length)
        .subtract(right.scale(f_w / 2))
        .subtract(up.scale(f_h / 2));

    return .{
        .w = w,
        .h = h,
        .origin = origin,
        .lower_left_corner = lower_left_corner,
        .right = right,
        .up = up
    };
}




// =====================================================================================
// =====================================================================================
// =====================================================================================
// =====================================================================================

fn loadColorTexture(allocator: std.mem.Allocator, gltf: Gltf, texture_info: ?Gltf.TextureInfo, _factor: []const f32) !stage3.Texture {
    const factor = vec3(_factor[0], _factor[1], _factor[2]);
    if (texture_info) |info| {
        const texture = gltf.data.textures.items[info.index];
        const image = gltf.data.images.items[texture.source.?];
        const im: *const zigimg.Image = @alignCast(@ptrCast(image.data.?));
        const data = try allocator.alloc(Vec3, im.width*im.height);
        var iter = im.iterator();
        while (iter.next()) |pixel| {
            data[iter.current_index-1] = vec3(pixel.r, pixel.g, pixel.b).mul(factor);
        }
        return .{
            .data = data,
            .w = @floatFromInt(im.width),
            .h = @floatFromInt(im.height),
            .w_int = im.width,
        };
    } else {
        const data = try allocator.alloc(Vec3, 1);
        data[0] = factor;
        return .{
            .data = data,
            .w = 1,
            .h = 1,
            .w_int = 1,
        };
    }
}

fn loadMaterial(allocator: std.mem.Allocator, gltf: Gltf, material_idx: Gltf.Index) !stage3.Material {
    const material = gltf.data.materials.items[material_idx];
    return .{
        .base_color = try loadColorTexture(allocator, gltf,
            material.metallic_roughness.base_color_texture,
            &material.metallic_roughness.base_color_factor),
        .emissive = try loadColorTexture(allocator, gltf,
            material.emissive_texture,
            &material.emissive_factor), // TODO: material.emissive_strength
    };
}

pub fn loadMaterials(gltf: Gltf, scene: *stage3.Scene) !void {
    const allocator = scene.arena.allocator();

    const materials = try allocator.alloc(stage3.Material, gltf.data.materials.items.len);
    for (materials, 0..) |*material, i| {
        material.* = try loadMaterial(allocator, gltf, i);
    }

    std.log.info("Materials count: {}", .{materials.len});

    scene.materials = materials;
}
