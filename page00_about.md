---
layout: page
title: About
description: Home
img: seg.png
permalink: index.html
sidebar: true
---

---

# {{site.data.about.title}}
By {{site.data.about.authors}}

{% for entry in site.data.about %}

{% if entry[0] != 'title' %}
{% if entry[0] != 'authors' %}
## {{entry[0]}}
{{entry[1]}}
{% endif %}
{% endif %}
{% endfor %}