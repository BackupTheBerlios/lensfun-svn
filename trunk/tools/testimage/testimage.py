#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""This program generates a bitmap which simulates lens errors.  Then, LensFun
can be used to rectify them.  This way, one can check whether LensFun is
working properly.  The test image contains vignetting and a rectangular grid.
The grid lines exhibit distortion and TCA.  The lens projection can be
non-rectilinear.

External dependencies: Python version 3, ImageMagick, exiv2.

The program takes the lens parameters from LensFun's database.  To the usual
search path it adds the current path (the xml files of which would override
those found in the usual locations).  It assumes an apect ratio of 3:2 for the
lens calibration.

As for the used camera, it takes the crop factor from the database.  For the
aspect ratio, it also assumes 3:2.

It is sensible to read the help output of this program (the --help option).

The following things are particularly interesting to check:

1. Vignetting, TCA, and distortion need to be corrected fully (apart from
   rounding errors)

    Note that for seeing the perfect grid, you must transform fisheyes into
    rectilinear.

2. Vignetting has to be corrected first, then TCA, then distortion.

    This is because this is the way they are calibrated: Vignetting is
    calibrated on a totally uncorrected image.  TCA is so, too, however, TCA
    doesn't care about vignetting noticably.  Well, and then the image is
    un-distorted.

    Note that LensFun seems to work the other way round.  But this is only
    because it uses the corrected destination image as the starting point: It
    uses the distortion formula to find the pixel on the distorted sensor
    image, then it applies the TCA correction formula to find the subpixel
    position, and finally, it applies the vignetting.

3. If he crop factor of the camera does not match the one of the calibration,
   it must work nevertheless.

    This is a tricky one.  If you have a sensor smaller than the one used for
    calibration, this program testimage.py has an easy job: It just has to use
    the corresponsing inner part of the resulting picture.

    However, for LensFun, it is more difficult because it works on the
    destination image.  Thus, cropping this image by the crop factor may lead
    to a wrong scaling of the distortion field relatively to the image.

    Let the distortion function be D(r).

    There are two coordinate systems, one for "r_u" in the undistorted
    destination image, and one for r_d = D(r_u) in the distorted sensor image.
    LensFun 0.2.8 scales the r_u by the ratio of the crop factors R_cf < 1,
    i.e. r_u' = r_u · R_cf.  But it scales r_d (i.e., the sensor image) the
    same way and this is wrong because it must use r_d' = r_d · D(R_cf).

    BTW, the scaling of r_u is arbitrary.  It just scales the resulting
    picture.  The above choice of R_cf means that one is probably close to the
    frame size so that no or little autoscaling is necessary.

    It is important to check this also for models with D(1) ≠ 1, e.g. poly5.

4. If the aspect ratio of the camera does not match the one of the calibration,
   it must work nevertheless.

    Classical case: Calibration with APS-C, but to-be-corrected image taken
    with a Four-Thirds sensor.


Todos:

- Support for non-rectilinear lens types
"""

import array, subprocess, math, os, argparse, sys
from xml.etree import ElementTree


parser = argparse.ArgumentParser(description="Generate test images for LensFun.")
parser.add_argument("lens_model_name", metavar="lens", help="Lens model name.  Must match an entry in LensFun's database exactly.")
parser.add_argument("camera_model_name", metavar="camera",
                    help="Camera model name.  Must match an entry in LensFun's database exactly.")
parser.add_argument("focal_length", metavar="focal length", type=float,
                    help="Focal length in mm.  Must match an entry in LensFun's database exactly, no interpolation is done.")
parser.add_argument("aperture", type=float,
                    help="Aperture in f-stops.  Must match an entry in LensFun's database exactly, no interpolation is done.")
parser.add_argument("distance", type=float,
                    help="Distance in metres.  Must match an entry in LensFun's database exactly, no interpolation is done.")
parser.add_argument("--width", type=int, default=600, help="The resulting bitmap's long edge width in pixels.  Default: 600")
parser.add_argument("--aspect-ratio", type=float, default=1.5, help="The resulting bitmap's aspect ratio, >= 1.  Default: 1.5")
parser.add_argument("--portrait", action="store_true",
                    help="Whether the resulting bitmap should be in portrait orientation.  Default: landscape")
parser.add_argument("--outfile", default="testimage.tiff", help="Path to the output file.  Default: testimage.tiff")
parser.add_argument("--no-vignetting", dest="vignetting", action="store_false",
                    help="Supresses simulation of vignetting.  *Much* faster.")
args = parser.parse_args()
lens_model_name, camera_model_name, focal_length, aperture, distance = \
                                args.lens_model_name, args.camera_model_name, args.focal_length, args.aperture, args.distance
focal_length, aperture, distance = float(focal_length), float(aperture), float(distance)
width, aspect_ratio, portrait = args.width, args.aspect_ratio, args.portrait


def get_database_elements():
    lens_element = camera_element = None
    distortion_element = vignetting_element = tca_element = None
    def crawl_directory(dirpath):
        nonlocal lens_element, camera_element, distortion_element, vignetting_element, tca_element
        for root, __, filenames in os.walk(dirpath):
            for filename in filenames:
                if filename.endswith(".xml"):
                    tree = ElementTree.parse(os.path.join(root, filename)).getroot()
                    for element in tree:
                        if camera_element is None and \
                           element.tag == "camera" and element.find("model").text == camera_model_name:
                            camera_element = element
                        elif element.tag == "lens" and element.find("model").text == lens_model_name:
                            if lens_element is not None:
                                print("This program cannot handle lens model names that occur multiple times in the\n"
                                      "read XML files.  Please use another lens.")
                                sys.exit(1)
                            lens_element = element
                            for calibration_element in element.find("calibration"):
                                if calibration_element.tag == "distortion" and \
                                   float(calibration_element.attrib["focal"]) == focal_length:
                                    distortion_element = calibration_element
                                elif calibration_element.tag == "tca" and \
                                   float(calibration_element.attrib["focal"]) == focal_length:
                                    tca_element = calibration_element
                                elif calibration_element.tag == "vignetting" and \
                                   float(calibration_element.attrib["focal"]) == focal_length and \
                                   float(calibration_element.attrib["aperture"]) == aperture and \
                                   float(calibration_element.attrib["distance"]) == distance:
                                    vignetting_element = calibration_element
    for path in [".", os.path.expanduser("~/.local/share/lensfun"), "/usr/share/lensfun", "/usr/local/share/lensfun"]:
        crawl_directory(path)
    if lens_element is None:
        print("Lens model name not found.")
        sys.exit(1)
    if camera_element is None:
        print("Camera model name not found.")
        sys.exit(1)
    return lens_element, camera_element, distortion_element, vignetting_element, tca_element

lens_element, camera_element, distortion_element, vignetting_element, tca_element = get_database_elements()


camera_cropfactor = float(camera_element.find("cropfactor").text)
lens_cropfactor = float(lens_element.find("cropfactor").text)
R_cf = lens_cropfactor / camera_cropfactor

def get_lens_aspect_ratio(lens_element):
    try:
        aspect_ratio = lens_element.find("aspect-ratio").text
    except AttributeError:
        return 1.5
    if ":" in aspect_ratio:
        numerator, denominator = aspect_ratio.split(":")
        return int(numerator) / int(denominator)
    return float(aspect_ratio)
lens_aspect_ratio = get_lens_aspect_ratio(lens_element)


def get_float_attribute(element, attribute_name, default=0):
    try:
        return float(element.attrib[attribute_name])
    except KeyError:
        return default

full_frame_diagonal = math.sqrt(36**2 + 24**2)

def get_projection_function():
    scaling = full_frame_diagonal / camera_cropfactor / focal_length
    def projection(r_vignetting):
        return math.atan(r_vignetting * scaling)
    lens_type = lens_element.get("type", "rectilinear")
    if lens_type == "stereographic":
        def projection(r_vignetting):
            return 2 * math.atan(r_vignetting * scaling / 2)
    elif lens_type == "fisheye":
        def projection(r_vignetting):
            return r_vignetting * scaling
    elif lens_type == "panoramic":
        assert False, "Panoramic lenses are not yet supported."
    elif lens_type == "equirectangular":
        assert False, "Equirectangular lenses are not yet supported."
    elif lens_type == "orthographic":
        def projection(r_vignetting):
            return math.asin(r_vignetting * scaling)
    elif lens_type == "equisolid":
        def projection(r_vignetting):
            return 2 * math.asin(r_vignetting * scaling / 2)
    elif lens_type == "fisheye_thoby":
        assert False, "Thoby fisheye lenses are not yet supported."
    return projection
projection = get_projection_function()
    
def get_distortion_function():
    def distortion(r):
        return 1
    if distortion_element is not None:
        model = distortion_element.attrib["model"]
        if model == "ptlens":
            a = get_float_attribute(distortion_element, "a")
            b = get_float_attribute(distortion_element, "b")
            c = get_float_attribute(distortion_element, "c")
            def distortion(r):
                return a * r**4 + b * r**3 + c * r**2 + (1 - a - b - c) * r
        elif model == "poly3":
            k1 = get_float_attribute(distortion_element, "k1")
            def distortion(r):
                return k1 * r**3 + (1 - k1) * r
        elif model == "fov1":
            ω = get_float_attribute(distortion_element, "omega")
            def distortion(r):
                return math.tan(r * ω) / (2 * tan(ω / 2))
        elif model == "poly5":
            k1 = get_float_attribute(distortion_element, "k1")
            k2 = get_float_attribute(distortion_element, "k2")
            def distortion(r):
                return r * (1 + k1 * r**2 + k2 * r**4)
    return distortion
distortion = get_distortion_function()

def get_tca_functions():
    def tca_red(r):
        return 1
    def tca_blue(r):
        return 1
    if tca_element is not None:
        model = tca_element.attrib["model"]
        if model == "linear":
            kr = get_float_attribute(tca_element, "kr", 1)
            def tca_red(r):
                return r * kr
            kb = get_float_attribute(tca_element, "kb", 1)
            def tca_blue(r):
                return r * kb
        elif model == "poly3":
            br = get_float_attribute(tca_element, "br")
            cr = get_float_attribute(tca_element, "cr")
            vr = get_float_attribute(tca_element, "vr", 1)
            def tca_red(r):
                return r**3 * br + r**2 * cr + r * vr
            bb = get_float_attribute(tca_element, "bb")
            cb = get_float_attribute(tca_element, "cb")
            vb = get_float_attribute(tca_element, "vb", 1)
            def tca_blue(r):
                return r**3 * bb + r**2 * cb + r * vb
    return tca_red, tca_blue
tca_red, tca_blue = get_tca_functions()

def get_vignetting_function():
    def vignetting(r):
        return 1
    if vignetting_element is not None:
        model = vignetting_element.attrib["model"]
        if model == "pa":
            k1 = get_float_attribute(vignetting_element, "k1")
            k2 = get_float_attribute(vignetting_element, "k2")
            k3 = get_float_attribute(vignetting_element, "k3")
            def vignetting(r):
                return 1 + k1 * r**2 + k2 * r**4 + k3 * r**6
    return vignetting
vignetting = get_vignetting_function()


class Image:

    def __init__(self, width, height):
        """width and height in pixels."""
        self.pixels = array.array("H", width * height * 3 * [16383])
        self.width = width
        self.height = height
        self._half_height = height / 2
        self._half_width = width / 2
        self._half_diagonal = math.sqrt(width**2 + height**2) / 2
        aspect_ratio = width / height
        aspect_ratio_correction = math.sqrt((lens_aspect_ratio**2 + 1) / (aspect_ratio**2 + 1))
        self.ar_plus_cf_correction = aspect_ratio_correction * R_cf

    def add_to_pixel(self, x, y, weighting, red, green, blue):
        """Adds the given brightness (from 0 to 1) to the pixel at the position x, y.
        `weighting` is multiplied with `red`, `green`, and `blue` before it is
        added, so it is some sort of opacity.  Coordinates must be integers.
        Coordinates outside the picture are ignored.
        """
        if 0 < x < self.width and 0 < y < self.height:
            offset = 3 * (y * width + x)
            weighting *= 65535
            self.pixels[offset] = max(min(self.pixels[offset] + int(red * weighting), 65535), 0)
            self.pixels[offset + 1] = max(min(self.pixels[offset + 1] + int(green * weighting), 65535), 0)
            self.pixels[offset + 2] = max(min(self.pixels[offset + 2] + int(blue * weighting), 65535), 0)

    def add_to_position(self, x, y, red, green, blue):
        """Adds the given brightness (from 0 to 1) to the pixel at the position x, y.
        The special bit here is that x and y are floats.
        """
        floor_x = int(math.floor(x))
        floor_y = int(math.floor(y))
        ceil_x = int(math.ceil(x)) if x != floor_x else int(x) + 1
        ceil_y = int(math.ceil(y)) if y != floor_y else int(y) + 1
        self.add_to_pixel(floor_x, floor_y, (ceil_x - x) * (ceil_y - y), red, green, blue)
        self.add_to_pixel(floor_x, ceil_y, (ceil_x - x) * (y - floor_y), red, green, blue)
        self.add_to_pixel(ceil_x, floor_y, (x - floor_x) * (ceil_y - y), red, green, blue)
        self.add_to_pixel(ceil_x, ceil_y, (x - floor_x) * (y - floor_y), red, green, blue)

    def r_vignetting(self, x, y):
        """Returns the r coordinate to x, y, for vignetting, i.e. r = 1 is
        half-diagonal.
        """
        x = (x - self._half_width) / self._half_diagonal
        y = (y - self._half_height) / self._half_diagonal
        return math.sqrt(x**2 + y**2) * R_cf

    def set_vignetting(self, function):
        for y in range(self.height):
            for x in range(self.width):
                offset = 3 * (y * self.width + x)
                for index in range(3):
                    self.pixels[offset + index] = \
                            max(min(int(self.pixels[offset + index] * function(self.r_vignetting(x, y))), 65535), 0)

    def rotate_by_90_degrees(self):
        """This changes the orientation to portrait.  Use this method shortly before
        writing the image, because further operations (creating grid etc) are
        undefined.
        """
        new_pixels = array.array("H", self.width * self.height * 3 * [0])
        for y in range(self.height):
            for x in range(self.width):
                offset_source = 3 * (y * self.width + x)
                offset_destination = 3 * ((self.width - x - 1) * self.height + y)
                new_pixels[offset_destination] = self.pixels[offset_source]
                new_pixels[offset_destination + 1] = self.pixels[offset_source + 1]
                new_pixels[offset_destination + 2] = self.pixels[offset_source + 2]
        self.pixels = new_pixels
        self.width, self.height = self.height, self.width

    def write(self, filepath):
        """Writes the image to a file.  If it is a TIFF or a JPEG, EXIF information is
        written, too.  ImageMagick and exiv2 are needed for anything beyond PPM
        output.
        """
        self.pixels.byteswap()
        ppm_data = "P6\n{} {}\n65535\n".format(self.width, self.height).encode("ascii") + self.pixels.tostring()
        self.pixels.byteswap()
        if filepath.endswith(".ppm"):
            open(filepath, "wb").write(ppm_data)
        else:
            subprocess.Popen(
                ["convert", "ppm:-", "-set", "colorspace", "RGB", filepath], stdin=subprocess.PIPE).communicate(ppm_data)
        if os.path.splitext(filepath)[1].lower() in [".tiff", ".tif", ".jpeg", "jpg"]:
            subprocess.call(["exiv2",
                             "-Mset Exif.Photo.FocalLength {}/10".format(int(focal_length * 10)),
                             "-Mset Exif.Photo.FNumber {}/10".format(int(aperture * 10)),
                             "-Mset Exif.Photo.SubjectDistance {}/10".format(int(distance * 10)),
                             "-Mset Exif.Photo.LensModel {}".format(lens_model_name),
                             "-Mset Exif.Image.Model {}".format(camera_model_name),
                             filepath])

    def create_grid(self, distortion, tca_red, tca_blue):
        line_brightness = 0.4
        def set_pixel(x, y):
            """y goes from -1 to 1, x domain depends on aspect ratio."""
            r = math.sqrt(x**2 + y**2) * self.ar_plus_cf_correction
            if r == 0:
                self.add_to_position(self._half_width, self._half_height, 1, 1, 1)
            else:
                r_distorted = distortion(r)
                scaling = (r_distorted / r) * self._half_height / R_cf
                self.add_to_position(x * scaling + self._half_width, y * scaling + self._half_height, 0, line_brightness, 0)
                scaling = tca_red(r_distorted) / r_distorted * (r_distorted / r) * self._half_height / R_cf
                self.add_to_position(x * scaling + self._half_width, y * scaling + self._half_height, line_brightness, 0, 0)
                scaling = tca_blue(r_distorted) / r_distorted * (r_distorted / r) * self._half_height / R_cf
                self.add_to_position(x * scaling + self._half_width, y * scaling + self._half_height, 0, 0, line_brightness)
        number_of_lines = 30
        for i in range(number_of_lines):
            points_per_line = self.width
            y = i * (2 / number_of_lines) - 1
            for j in range(points_per_line):
                x = j * (2 * lens_aspect_ratio / points_per_line) - lens_aspect_ratio
                set_pixel(x, y)
            points_per_line = self.height
            x = i * (2 * lens_aspect_ratio / number_of_lines) - lens_aspect_ratio
            for j in range(points_per_line):
                y = j * (2 / points_per_line) - 1
                set_pixel(x, y)


image = Image(width, int(width / aspect_ratio))
image.create_grid(distortion, tca_red, tca_blue)
if args.vignetting:
    image.set_vignetting(vignetting)
if portrait:
    image.rotate_by_90_degrees()
image.write(args.outfile)
