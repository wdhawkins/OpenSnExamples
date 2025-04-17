seed = 48; // Change for a different city layout
block_size = 100;
num_blocks = 2;
street_width = 10;
sidewalk_width = 5;
building_min_size = 10;
building_max_size = 40;
height_min = 20;
height_max = 100;
overall_size = num_blocks * (block_size + street_width) + sidewalk_width * 2;
ground_height = 1;

module rectangular_building(width, depth, height) {
    color("gray") translate([0, 0, ground_height]) cube([width, depth, height]);
}

module city_block(x_offset, y_offset) {
    lot_size = block_size / 3;
    for (i = [0:2]) {
        for (j = [0:2]) {
            x = x_offset + i * lot_size;
            y = y_offset + j * lot_size;
            width = rands(building_min_size, building_max_size, 1, seed + x)[0];
            depth = rands(building_min_size, building_max_size, 1, seed + y)[0];
            height = rands(height_min, height_max, 1, seed + x * y)[0];
            
            translate([x + (lot_size - width) / 2, y + (lot_size - depth) / 2, 0]) {
                rectangular_building(width, depth, height);
            }
        }
    }
}

module city_layout() {
    color("lightgray") translate([-sidewalk_width, -sidewalk_width, 0]) cube([overall_size, overall_size, ground_height]);
    for (i = [0:num_blocks-1]) {
        for (j = [0:num_blocks-1]) {
            x_offset = i * (block_size + street_width) + sidewalk_width;
            y_offset = j * (block_size + street_width) + sidewalk_width;
            city_block(x_offset, y_offset);
        }
    }
}

module source() {
    x = rands(sidewalk_width, overall_size - sidewalk_width, 1, seed + 300)[0];
    y = rands(sidewalk_width, overall_size - sidewalk_width, 1, seed + 400)[0];
    color("red") translate([x, y, ground_height]) cube([1, 1, 1]);
}

module air_cube() {
    color("cyan", 0.3) translate([-sidewalk_width, -sidewalk_width, 0]) 
        cube([overall_size, overall_size, height_max]);
}

union() {
    city_layout();
}

source();

air_cube();