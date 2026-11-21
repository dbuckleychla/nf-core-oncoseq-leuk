"""Render portable HTML with compact images and complete selected evidence."""
import base64
import hashlib
import io
import json
import re
import warnings
import xml.etree.ElementTree as ET
from pathlib import Path

from jinja2 import Environment, FileSystemLoader, select_autoescape
from PIL import Image

from adapters import chromosome_key

Image.MAX_IMAGE_PIXELS = 160_000_000


def image_uri(path):
    path = Path(path)
    if path.suffix.lower() == '.svg':
        tree = ET.fromstring(path.read_bytes())
        for parent in tree.iter():
            for child in list(parent):
                if child.tag.split('}')[-1] in ('script', 'foreignObject'):
                    parent.remove(child)
            for key in list(parent.attrib):
                if key.lower().startswith('on') or (key.split('}')[-1] == 'href'
                                                     and not parent.attrib[key].startswith('#')):
                    del parent.attrib[key]
        raw = ET.tostring(tree)
        mime = 'image/svg+xml'
    else:
        with Image.open(path) as image:
            image.thumbnail((3000, 5000), Image.Resampling.LANCZOS)
            buffer = io.BytesIO()
            if image.mode=='RGBA' or 'transparency' in image.info:
                rgba=image.convert('RGBA')
                rgb=Image.new('RGB',rgba.size,'white')
                rgb.paste(rgba,mask=rgba.getchannel('A'))
            else:rgb=image.convert('RGB')
            rgb.save(buffer, format='WEBP', lossless=True, method=6)
            raw, mime = buffer.getvalue(), 'image/webp'
    return 'data:' + mime + ';base64,' + base64.b64encode(raw).decode()


def render(input, output, template, figures=None):
    source = Path(input)
    data = json.loads(source.read_text())
    if data.get('schema_version') != 1 or not data.get('sample'):
        raise ValueError('Unsupported or incomplete review data')
    assets = list(data.get('figures', []))
    if figures:
        root = Path(figures)
        manifest = root / 'figures.json'
        if manifest.exists():
            assets.extend(json.loads(manifest.read_text()))
            for item in assets:
                if item.get('generated'):
                    item['path'] = str(root / item['path'])
    encoded, gallery = {}, []
    for item in assets:
        path = Path(item['path'])
        if not path.is_absolute():
            path = source.parent / path
        if not path.is_file():
            raise ValueError('Missing report figure: ' + str(path))
        key = hashlib.sha256(path.read_bytes()).hexdigest()
        if key not in encoded:
            encoded[key] = image_uri(path)
        match = re.search(r'(?:^|_)chr([0-9]+|X|Y|M)(?:_|\.)', item.get('name', path.name))
        order = chromosome_key(match[1]) if match else (27, item.get('name', path.name))
        gallery.append(dict(item, image_key=key, order=order))
    gallery.sort(key=lambda r: (r['category'], r['order'], r['name']))
    data['figures'] = gallery
    data['images'] = encoded
    payload = json.dumps(data, separators=(',', ':'), allow_nan=False).replace('<', '\\u003c').replace('&', '\\u0026')
    environment = Environment(loader=FileSystemLoader(str(Path(template).parent)), autoescape=select_autoescape(['html']))
    html = environment.get_template(Path(template).name).render(sample=data['sample'], data_json=payload)
    target = Path(output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(html)
    if target.stat().st_size > 50 * 1024 * 1024:
        warnings.warn('Review HTML exceeds 50 MiB; all selected evidence was retained')
    return target
