// Parameters
main_text = "OPENS";    // Main text
subscript_text = "n";   // Subscript
font_size_main = 20;    // Main text font size
font_size_sub = 12;     // Subscript font size
extrude_depth = 10;     // Thickness of the letters
font_style = "Liberation Sans:style=Bold";  // Font selection

// Generate the main text (Extruded)
module main_text_model() {
    linear_extrude(height=extrude_depth)
        text(main_text, size=font_size_main, font=font_style, halign="center", valign="center");
}

// Generate the subscript text (Extruded)
module subscript_model() {
    translate([font_size_main * 1.8 + 17, -font_size_main * 0.25 - 5, 0])  // Adjust subscript position
    linear_extrude(height=extrude_depth)
        text(subscript_text, size=font_size_sub, font=font_style, halign="center", valign="center");
}

module block() {
    color("cyan", 0.3) translate([-52, -22, 0]) cube([115, 40, 10]);
}

// Render the fully extruded model
main_text_model();
subscript_model();
block();