---
layout: page
title: Code and Data
img: seg.png 
permalink: code
sidebar: true
---

---

{% if site.data.code %}
## Code
{% for script in site.data.code %}
* [**{{script.name}}**]({{site.url}}/{{site.baseurl}}/software/{{script.name}})
  \| {{script.desc}}
{% endfor %}
{% endif %}

{% if site.data.datasets %}
## Data Sets
{% for ds in site.data.datasets %}
* [{{ds.name}}]({%if ds.storage !=
  'remote'%}{{site.url}}/{{site.baseurl}}/datasets/{{ds.link}}{%
  else%}{{site.link}}{% endif %}) \| {% if ds.filetype %}(filetype:
  {{ds.filetype}}){%endif%}{% if ds.filesize %}({{ds.filesize}}){%endif%}{%
  if ds.storage == 'remote' %} DOI: {{ds.DOI}}{%endif%}
{% endfor %}
{% endif %}

{% if site.data.figures %}
## Main Text Figure Generation

{% for fig in site.data.figures %}
<article class="post">

<a class="post-thumbnail" style="background-image: url({{site.url}}/{{site.baseurl}}/assets/img/{{fig.pic}})" href=""> </a>

<div class="post-content">
<b class="post-title"><a href="{{site.url}}/{{site.baseurl}}/software/{{fig.filename}}">{{fig.title}}</a></b>
<p> {{fig.desc}}</p>

<i>Necessary Data Sets </i><br/>
{% for ds in fig.req %}
{% if ds.storage == 'local' %}
{% assign link = "{{site.url}}/{{site.baseurl}}/datasets/{{ds.link}}" %}
{% else %}
{% assign link = "{{ds.link}}" %}
{% endif %}
• <a style="font-size: 0.9em;" href="{{link}}">  {{ds.title}} </a><br/>
{% endfor %}
</div>
</article>
{%endfor%}
{% endif %}



## SI Text Figure Generation
{% for fig in site.data.si_figures %}
<article class="post">

<a class="post-thumbnail" style="background-image: url({{site.url}}/{{site.baseurl}}/assets/img/{{fig.pic}})" href="{{site.baseurl}}/figures/{{fig.pdf}}"> </a>

<div class="post-content">
<b class="post-title"><a href="{{site.url}}/{{site.baseurl}}/software/{{fig.filename}}">{{fig.title}}</a></b>
<p> {{fig.desc}}</p>

<i>Necessary Data Sets </i><br/>
{% assign label = 0 %}
{% for ds in fig.req %}
{% if ds.storage == 'local' %}
{% assign link = "{{site.url}}/{{site.baseurl}}/datasets/{{ds.link}}" %}
{% else %}
{% assign link = "{{ds.link}}" %}
{% endif %}
{% if ds.stan == true %}
{% if label == 0 %}
<i> Necessary Stan Models </i><br/>
{% assign label = 1 %}
{% endif %}
• <a style="font-size: 0.9em;" href="{{site.url}}/{{site.baseurl}}/software/{{ds.link}}"> {{ds.title}}</a><br/>

{% else %}
• <a style="font-size: 0.9em;" href="{{link}}"> {{ds.title}} </a> <br />
{% endif %}
{% endfor %}

</div>
</article>
{%endfor%}
