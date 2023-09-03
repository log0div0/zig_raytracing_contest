
const std = @import("std");

const main = @import("main.zig");
const linalg = @import("linalg.zig");
const stage3 = @import("stage3.zig");

const Vec3 = linalg.Vec3;
const Grid = linalg.Grid;
const Bbox = linalg.Bbox;

// =====================================================================================
// =====================================================================================
// =====================================================================================
// =====================================================================================

pub const Triangle = struct {
	centroid: [3]f32,
	bbox: Bbox,
	pos: [3]Vec3,
	data: stage3.Triangle.Data,
};


const BVHNode = struct {
	bbox: Bbox,
	payload: union(enum) {
		leaf: struct {
			begin: usize, // TODO: fix type
			end: usize, // TODO: fix type
		},
		interior: struct {
			second_child_offset: usize, // TODO: fix type
		},
	},
};

const Cell = struct {
	triangles: []usize = &.{},
	num_triangles: usize = 0,
    bvh_nodes: []BVHNode = &.{},
    num_bvh_nodes: usize = 0,
};

pub const Geometry = struct {
	triangles: std.MultiArrayList(Triangle).Slice,
	grid: Grid,
	cells: []Cell,
	arena: std.heap.ArenaAllocator,

	pub fn init(allocator: std.mem.Allocator) Geometry {
	    var geometry: Geometry = undefined;
	    geometry.arena = std.heap.ArenaAllocator.init(allocator);
	    return geometry;
	}

	pub fn deinit(self: *Geometry) void {
	    self.arena.deinit();
	}

	fn initGrid(geometry: *Geometry) !void {
	    var mean_triangle_size = Vec3.zeroes();
	    var bbox: Bbox = .{};

	    for (geometry.triangles.items(.bbox)) |triangle_bbox| {
            mean_triangle_size = Vec3.add(mean_triangle_size, triangle_bbox.size());
            bbox = bbox.unionWith(triangle_bbox);
	    }

	    mean_triangle_size = mean_triangle_size.div(Vec3.fromScalar(@floatFromInt(geometry.triangles.len)));
	    const mean_triangle_count = bbox.size().div(mean_triangle_size);

	    std.log.info("Mean triangle count: {d:.1}", .{mean_triangle_count.data});
	    const resolution = main.config.grid_resolution orelse mean_triangle_count.div(Vec3.fromScalar(4)).ceil().toInt(u32);
	    std.log.info("Grid resolution: {}", .{resolution.data});

	    geometry.grid = Grid.init(bbox, resolution);
	}

	fn initCells(geometry: *Geometry) !void {
		const grid = &geometry.grid;
		const num_cells = grid.resolution.reduceMul();
		const allocator = geometry.arena.allocator();
		const cells = try allocator.alloc(Cell, num_cells);
	    @memset(cells, .{});
	    for (geometry.triangles.items(.pos)) |pos| {
            const min = grid.getCellPos(Vec3.min(pos[0], Vec3.min(pos[1], pos[2])));
            const max = grid.getCellPos(Vec3.max(pos[0], Vec3.max(pos[1], pos[2])));

            for (min.z()..max.z()+1) |z| {
                for (min.y()..max.y()+1) |y| {
                    for (min.x()..max.x()+1) |x| {
                        const index = grid.getCellIdx(x,y,z);
                        cells[index].num_triangles += 1;
                    }
                }
            }
	    }
	    geometry.cells = cells;
	    var min_triangles: usize = std.math.maxInt(usize);
	    var max_triangles: usize = 0;
	    var empty_cells: usize = 0;
	    var total_triangles_count: usize = 0;
	    for (cells) |*cell| {
	        if (cell.num_triangles != 0) {
	            min_triangles = @min(min_triangles, cell.num_triangles);
	            max_triangles = @max(max_triangles, cell.num_triangles);
	        } else {
	            empty_cells += 1;
	        }
	        cell.triangles = try allocator.alloc(usize, cell.num_triangles);
	        total_triangles_count += cell.num_triangles;
	        cell.num_triangles = 0;
	    }
	    const mean_triangles = total_triangles_count / (grid.numCells() - empty_cells);
	    std.log.info("Empty cells: {}/{} ({d:.2}%) min triangles: {} max triangles: {} mean_triangles: {}",
	        .{empty_cells, grid.numCells(),
	            @as(f32, @floatFromInt(empty_cells)) / @as(f32, @floatFromInt(grid.numCells())) * 100,
	            min_triangles, max_triangles, mean_triangles});

	    for (geometry.triangles.items(.pos), 0..) |pos, triangle_index| {
            const min = grid.getCellPos(Vec3.min(pos[0], Vec3.min(pos[1], pos[2])));
            const max = grid.getCellPos(Vec3.max(pos[0], Vec3.max(pos[1], pos[2])));

            for (min.z()..max.z()+1) |z| {
                for (min.y()..max.y()+1) |y| {
                    for (min.x()..max.x()+1) |x| {
                        const cell_index = grid.getCellIdx(x,y,z);
                        var cell = &geometry.cells[cell_index];
                        cell.triangles[cell.num_triangles] = triangle_index;
                        cell.num_triangles += 1;
                    }
                }
            }
	    }

	    std.log.info("Unique triangle count: {}/{} ({d:.2}%)",
	        .{geometry.triangles.len, total_triangles_count,
	            @as(f32, @floatFromInt(geometry.triangles.len)) / @as(f32, @floatFromInt(total_triangles_count)) * 100});
	}

	const SortCtx = struct {
		dim: u8,
		centroids: [][3]f32,
		indices: []usize,

		pub fn lessThan(ctx: @This(), a: usize, b: usize) bool {
		    return ctx.centroids[ctx.indices[a]][ctx.dim] < ctx.centroids[ctx.indices[b]][ctx.dim];
		}

		pub fn swap(ctx: @This(), a: usize, b: usize) void {
		    return std.mem.swap(usize, &ctx.indices[a], &ctx.indices[b]);
		}
	};

	fn recursiveBuildBVH(geometry: Geometry, cell: *Cell, begin: usize, end: usize) void
	{
		const num_triangles = end - begin;
		const current_node_idx = cell.num_bvh_nodes;
		cell.num_bvh_nodes += 1;

		var bbox: Bbox = .{};
		for (begin..end) |i| {
			bbox = bbox.unionWith(geometry.triangles.items(.bbox)[cell.triangles[i]]);
		}

		if (num_triangles <= main.config.max_triangles_per_bvh_node) {
			cell.bvh_nodes[current_node_idx] = .{
				.bbox = bbox,
				.payload = .{
					.leaf = .{
						.begin = begin,
						.end = end,
					},
				},
			};
			return;
		}

		const dim = bbox.size().maxDim();
		std.mem.sortContext(begin, end, SortCtx{
			.dim = dim,
			.centroids = geometry.triangles.items(.centroid),
			.indices = cell.triangles,
		});

		const mid = (begin + end) / 2;

		const left_child_idx = cell.num_bvh_nodes;
		recursiveBuildBVH(geometry, cell, begin, mid);
		const right_child_idx = cell.num_bvh_nodes;
		recursiveBuildBVH(geometry, cell, mid, end);

		const left_child = &cell.bvh_nodes[left_child_idx];
		const right_child = &cell.bvh_nodes[right_child_idx];

		cell.bvh_nodes[current_node_idx] = .{
			.bbox = left_child.bbox.unionWith(right_child.bbox),
			.payload = .{
				.interior = .{
					.second_child_offset = right_child_idx,
				}
			},
		};
	}

	fn initBVH(geometry: *Geometry) !void {
		const allocator = geometry.arena.allocator();
		for (geometry.cells) |*cell| {
			if (cell.num_triangles != 0) {
				cell.bvh_nodes = try allocator.alloc(BVHNode, 2*cell.num_triangles - 1);
			}
		}
		for (geometry.cells) |*cell| {
			if (cell.num_triangles != 0) {
				recursiveBuildBVH(geometry.*, cell, 0, cell.num_triangles);
			}
		}
	}

	pub fn build(self: *Geometry) !void {
	    try self.initGrid();
	    try self.initCells();
	    try self.initBVH();
	}

	pub fn bakeInto(geometry: Geometry, scene: *stage3.Scene) !void {
	    const allocator = scene.arena.allocator();

	    var bvh_counter: usize = 0;
	    var triangle_counter: usize = 0;
	    const cells = try allocator.alloc(stage3.Cell, geometry.cells.len);
	    for (geometry.cells, cells) |src, *dst| {
	    	dst.* = .{
	    		.bvh_offset = if (src.num_triangles != 0) bvh_counter else std.math.maxInt(usize),
	    	};
	    	bvh_counter += src.num_bvh_nodes;
	    	triangle_counter += src.num_triangles;
	    }

	    const bvh_nodes = try allocator.alloc(stage3.BVHNode, bvh_counter);
	    var triangles = std.MultiArrayList(stage3.Triangle){};
	    try triangles.resize(allocator, triangle_counter);

	    bvh_counter = 0;
	    triangle_counter = 0;
	    for (geometry.cells, cells) |src_cell, dst_cell| {
	    	for (src_cell.bvh_nodes[0..src_cell.num_bvh_nodes]) |src_node| {
	    		const dst_node = &bvh_nodes[bvh_counter];
	    		switch (src_node.payload) {
	    			.leaf => |leaf| {
				    	dst_node.* = .{
				    		.bbox = src_node.bbox,
				    		.num_triangles = leaf.end - leaf.begin,
				    		.offset = .{
				    			.first_triangle_offset = triangle_counter + leaf.begin,
				    		}
				    	};
	    			},
	    			.interior => |interior| {
				    	dst_node.* = .{
				    		.bbox = src_node.bbox,
				    		.num_triangles = 0,
				    		.offset = .{
				    			.second_child_offset = dst_cell.bvh_offset + interior.second_child_offset,
				    		}
				    	};
	    			}
	    		}
		    	bvh_counter += 1;
	    	}
	    	for (src_cell.triangles[0..src_cell.num_triangles]) |triangle_index| {
		    	const src = geometry.triangles.get(triangle_index);
		    	const dst = stage3.Triangle{
		    	    .pos = stage3.Triangle.Pos.init(src.pos[0], src.pos[1], src.pos[2]),
		    	    .data = src.data,
		    	};
		    	triangles.set(triangle_counter, dst);
		    	triangle_counter += 1;
	    	}
	    }

	    scene.grid = geometry.grid;
	    scene.cells = cells;
	    scene.bvh_nodes = bvh_nodes;
	    scene.triangles_pos = triangles.items(.pos);
	    scene.triangles_data = triangles.items(.data);
	}
};
