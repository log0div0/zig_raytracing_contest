const std = @import("std");
const builtin = @import("builtin");
const zigimg = @import("zigimg");
const zigargs = @import("zigargs");

const linalg = @import("linalg.zig");
const stage1 = @import("stage1.zig");
const stage2 = @import("stage2.zig");
const stage3 = @import("stage3.zig");

test {
    _ = @import("linalg.zig");
}

const Vec3u = linalg.Vec3u;

// =====================================================================================
// =====================================================================================
// =====================================================================================
// =====================================================================================

fn getDuration(start_time: i128) @TypeOf(std.fmt.fmtDuration(1)) {
    const end_time = std.time.nanoTimestamp();
    return std.fmt.fmtDuration(@intCast(end_time - start_time));
}

pub const std_options = struct {
    pub const log_level = if (builtin.mode == .Debug) .debug else .info;
};

pub const CmdlineArgs = struct {
    in: []const u8 = "input.gltf",
    out: []const u8 = "output.png",
    width: ?u16 = null,
    height: ?u16 = null,
};

pub fn loadFile(path: []const u8, allocator: std.mem.Allocator) ![]const u8 {
    const file = try std.fs.cwd().openFile(path, .{});
    defer file.close();

    const file_size = try file.getEndPos();

    const buf = try allocator.alloc(u8, file_size);
    errdefer allocator.free(buf);

    const bytes_read = try file.readAll(buf);
    std.debug.assert(bytes_read == file_size);

    return buf;
}

const Config = struct {
    grid_resolution: ?Vec3u,
    max_triangles_per_bvh_node: usize,
    num_threads: ?u8,
    num_samples: u16,
    max_bounce: u16,

    fn load(path: []const u8, allocator: std.mem.Allocator) !Config {
        const buf = try loadFile(path, allocator);
        defer allocator.free(buf);
        var parsed_str = try std.json.parseFromSlice(Config, allocator, buf, .{});
        defer parsed_str.deinit();
        return parsed_str.value;
    }
};

pub var config: Config = undefined;

pub fn main() !void {
    std.log.info("{}", .{std_options.log_level});

    const start_time = std.time.nanoTimestamp();

    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(gpa.deinit() == .ok);
    const allocator = gpa.allocator();

    const args = try zigargs.parseForCurrentProcess(CmdlineArgs, allocator, .print);
    defer args.deinit();

    config = try Config.load("config.json", allocator);
    std.log.info("Num samples: {}, max bounce {}", .{config.num_samples, config.max_bounce});

    const num_threads = config.num_threads orelse try std.Thread.getCpuCount();
    const threads = try allocator.alloc(std.Thread, num_threads);
    defer allocator.free(threads);
    std.log.info("Num threads: {}", .{num_threads});

    var camera: stage3.Camera = undefined;

    var scene = stage3.Scene.init(allocator);
    defer scene.deinit();
    {
        var geometry = stage2.Geometry.init(allocator);
        defer geometry.deinit();
        {
            const loading_time = std.time.nanoTimestamp();
            var gltf = try stage1.loadGltfFile(allocator, args.options.in, threads);
            defer gltf.deinit();
            std.log.info("Loaded in {}", .{getDuration(loading_time)});

            const preprocessing_time = std.time.nanoTimestamp();
            camera = try stage1.loadCamera(gltf, args.options.width, args.options.height);
            try stage1.loadMaterials(gltf, &scene);
            try stage1.loadGeometry(gltf, &geometry);
            std.log.info("Preprocessed in {}", .{getDuration(preprocessing_time)});
        }

        const compile_time = std.time.nanoTimestamp();
        try geometry.build();
        try geometry.bakeInto(&scene);
        std.log.info("Compiled in {}", .{getDuration(compile_time)});
    }

    var img = try zigimg.Image.create(allocator, camera.w, camera.h, .rgb24);
    defer img.deinit();

    const render_time = std.time.nanoTimestamp();
    try scene.render(threads, camera, img);
    std.log.info("Rendered in {}", .{getDuration(render_time)});

    const save_time = std.time.nanoTimestamp();
    try img.writeToFilePath(args.options.out, .{.png = .{}});
    std.log.info("Saved in {}", .{getDuration(save_time)});

    std.log.info("Done in {}", .{getDuration(start_time)});
}
